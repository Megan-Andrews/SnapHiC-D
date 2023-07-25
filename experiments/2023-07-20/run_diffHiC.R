rm(list = ls())
gc() # collect garbage

library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

snap_dir = "/home/maa160/SnapHiC-D/"

library('rhdf5')
library(pryr)
#library(purrr)
library(InteractionSet)
library(dplyr)
library(data.table)
library(Matrix)
library(sparsesvd)
library(umap)
library(csaw)
library(multiHiCcompare)
library(pbapply)



read_hic_cool <- function(cond1_cooldir, cond2_cooldir, chr){
  coolpaths = lapply(list(cond1_cooldir, cond2_cooldir), function(p){
    file.path(p,list.files(p))})
  groups = unlist(lapply(c(1:2), function(i){
    rep(paste0('cond', i), length(coolpaths[[i]]))}))
  coolpaths = unlist(coolpaths)
  all_pixels = pblapply(c(1:length(coolpaths)), function(coolpaths, id){print(id)
    pixels <- h5read(coolpaths[id], c('pixels'));H5close()
    pixels = data.frame(pixels)
    pixels$chr = chr
    pixels[c('chr', 'bin1_id', 'bin2_id', 'count')]
  }, coolpaths = coolpaths)
  all_pixels
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one hic dataframe including all cells 

rbind_hic_dfs <- function(hic_dfs){
  all_hic_df = lapply(c(1:length(hic_dfs)), function(id){
    df = hic_dfs[[id]]
    df = df[,c('chr','bin1_id','bin2_id','count')]
    df$cell = paste0('cell',id)
    df
  })
  all_hic_df = rbindlist(all_hic_df)
  all_hic_df
}


# input: one hic dataframe including chr, bin1_id, bin2_id, count, and cell columns
# output: one hic dataframe with a shape of original one including BandNorm counts

bandnorm <- function(hic_df){
  hic_df$diag = hic_df$bin2_id - hic_df$bin1_id
  hic_df = hic_df[hic_df$diag != 0, ]
  band_info <- hic_df %>% group_by(diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(diag) %>% summarise(depth = mean(band_depth))
  hic_df <- hic_df %>% left_join(alpha_j, by = c("diag")) %>% 
    left_join(band_info, by = c("diag", "cell")) %>% 
    mutate(BandNorm = count/band_depth * depth) %>% 
    select(-c(band_depth, depth, count))
  #colnames(hic_df) = c('bin1_id', 'bin2_id', 'cell', 'diag', 'count')
  hic_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one InteractionSet object 

make_iset <- function(dfs){
  dfs = lapply(c(1:length(dfs)), function(data_frames, id){
    df = data_frames[[id]]
    colnames(df) = c('chr', 'bin1_id', 'bin2_id', paste0('count',id))
    data.table(df, key = c('chr', 'bin1_id', 'bin2_id'))
  }, data_frames = dfs)
  merged_df =  Reduce(function(dt1,dt2){merge.data.table(dt1,dt2,all = TRUE)}, 
                      dfs)
  merged_df[is.na(merged_df)] = 0
  chrs = merged_df$chr
  source_gr = GRanges(seqnames = merged_df$chr, 
                      ranges = IRanges(merged_df$bin1_id, end = merged_df$bin1_id+1))
  target_gr = GRanges(seqnames = merged_df$chr,
                      ranges = IRanges(merged_df$bin2_id, end = merged_df$bin2_id+1))
  GI = GInteractions(source_gr, target_gr)
  iset = InteractionSet(as(merged_df[,4:ncol(merged_df)], 'matrix'), GI)
  elementMetadata(iset) = chrs
  iset
}


# input: one InteractionSet object
# output: one InteractionSet that has been filtered to 
# exclude filter regions, and only include TSS regions
filter_regions_iset <- function(isets, filter_regions_path, resolution){  
  filter_regions = fread(filter_regions_path)
  filter_regions = filter_regions[,c(1:2)]
  colnames(filter_regions) = c('chr','start')
  filter_regions$bin_id = as.integer(filter_regions$start/resolution)
  filter_regions_dt = as.data.table(filter_regions)
  filter_regions_dt = unique(filter_regions_dt, by = c('chr','bin_id'))
  filter_regions_dt$include = rep(FALSE, length(filter_regions_dt$chr))
  
  chrs <- elementMetadata(isets)$X
  binIds <- anchorIds(isets, type="both")
  binId_dt <- data.table(bin1_id = binIds$first, bin2_id=binIds$second, chrs=chrs) 

  # Perform the join and filtering operation
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin1_id" = "bin_id", "chrs" ="chr"))
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin2_id" = "bin_id", "chrs" ="chr"))
  binId_dt[is.na(include.x), include.x := TRUE]
  binId_dt[is.na(include.y), include.y := TRUE]
  
  keep_filter_regions = binId_dt$include.x & binId_dt$include.y
  isets <- isets[keep_filter_regions,]
  
  return(isets)
}


# input: one InteractionSet object
# output: one InteractionSet that has been filtered to 
# exclude filter regions, and only include TSS regions
TSS_filter_iset <- function(iset, TSS_regions_path, resolution){  
  # # only including TSS regions 
  TSS_regions = fread(TSS_regions_path)
  TSS_regions$bin_id = as.integer((TSS_regions$start+TSS_regions$end)/(2*resolution))
  TSS_regions_dt = as.data.table(TSS_regions)
  TSS_regions_dt$include = rep(TRUE, length(TSS_regions_dt$chr))
  TSS_regions_dt = unique(TSS_regions_dt[,c("bin_id","include","chr")])
  
  chrs <- elementMetadata(isets)$X
  binIds <- anchorIds(isets, type="both")
  binId_dt <- data.table(bin1_id = binIds$first, bin2_id=binIds$second, chrs=chrs) 
  
  # Perform the join and filtering operation
  binId_dt = left_join(binId_dt, TSS_regions_dt, by = c("bin1_id" = "bin_id", "chrs" ="chr"))
  binId_dt = left_join(binId_dt, TSS_regions_dt, by = c("bin2_id" = "bin_id", "chrs" ="chr"))
  binId_dt[is.na(include.x), include.x := FALSE]
  binId_dt[is.na(include.y), include.y := FALSE]
  
  keep_gene_transcript = binId_dt$include.x & binId_dt$include.y
  isets <- isets[keep_gene_transcript,]

  return(isets)
}

### TODO: coarsen_iset?
# input: one InteractionSet
# input: n_groups, the number pseudo_bulk groups per cell type
coarsen_iset <- function(isets, n_groups){
  n_samples = dim(assay(isets))[2]
  if ((n_samples/2) %% n_groups != 0){
    print("pseudo-bulk will not have the same number of samples")
    return(NULL)
  }
  
  group_size = (n_samples/2) / n_groups
  for (i in seq(1, n_samples - group_size + 1, by = group_size)) {
    new_col <- rowSums(assay(isets)[, i:(i + group_size - 1)])
    col_name <- paste0("coarse_", i, "-", i + group_size - 1)
    col_names = colnames(isets)
    col_names[ceiling(i/group_size)]=col_name
    colnames(isets)=col_names
    assay(isets)[, ceiling(i/group_size)]= new_col
  }
  
  # Remove the original columns
  isets = isets[, 1:(n_samples/group_size)]
  return(isets)
}



chrom_size_filepath <- paste0(snap_dir, 'ext/hg19.chrom.sizes')
filter_regions_path <- paste0(snap_dir, "ext/hg19_filter_regions.txt")
TSS_regions_path <- paste0(snap_dir, "/ext/hg19.refGene.transcript.TSS.061421.txt")

result_directory <- paste0("/project/compbio-lab/scHi-C/Lee2019/results/2023-07-25/")

run_diffHic <- function(typeA_dir, typeB_dir, output_file, filter_criteria, batch_size){
    ##### BEGIN DIFFHIC #####
    print("loading data")

    dfs = read_hic_cool(typeA_dir, typeB_dir, "chr22")
    isets = make_iset(dfs)
    isets = isets[,c(1:90, (dim(isets)[2]-90+1):dim(isets)[2])] # get first 90 (cond1) and last 90 (cond2) samples
    if(batch_size != 90){
        isets = coarsen_iset(isets, batch_size)
    }

    assayNames(isets) = 'counts'
    anchors = anchors(isets,type="both")
    GI = GInteractions(anchors$first, anchors$second)
    y <- asDGEList(isets)

    y <- normOffsets(y, se.out=TRUE) # normalization

    counts = y$counts
    rownames(counts) <- NULL
    colnames(counts) <- NULL
    chrs = elementMetadata(isets)$X
    isets = InteractionSet(as(counts, "matrix"), GI) # recreate InteractionSet
    elementMetadata(isets) = chrs
    assayNames(isets) = 'counts'

    print("filtering uninteresting bin pairs")  
    isets = filter_regions_iset(isets, filter_regions_path, 100000)
    
    if(filter_criteria == "A"){
        keep <- rowSums(assay(isets)>0) > ncol(assay(isets))*0.1 # filter A
    }else if (filter_criteria == "B"){
        keep <- aveLogCPM(asDGEList(isets)) > 2 # filter B
    }

    isets <- isets[keep,]
    y <- asDGEList(isets)

    # modeling and testing 
    group = factor(c(rep(1,batch_size) , rep(2, batch_size)))
    design = model.matrix(~group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    result <- glmQLFTest(fit)
    adj.p <- p.adjust(result$table$PValue, method="BH")
    sum(adj.p <= 0.05)
    useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
    inter.frame <- as.data.frame(interactions(isets))[,useful.cols]
    results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
    o.r <- order(results.r$PValue)
    results.r = results.r[o.r,]

    write.csv(results.r, file = paste0(output_file), row.names = FALSE)
}

typeA_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_21yr/cond1"
typeB_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_21yr/sim_cond2"

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_90x90_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_90x90_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 90)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 90)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_45x45_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_45x45_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 45)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 45)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_30x30_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_30x30_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 30)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 30)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_18x18_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_18x18_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 18)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 18)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_15x15_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_15x15_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 15)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 15)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_10x10_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_10x10_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 10)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 10)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_9x9_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_9x9_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 9)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 9)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_6x6_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_6x6_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 6)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 6)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_5x5_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_5x5_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 5)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 5)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_3x3_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_3x3_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 3)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 3)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_2x2_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_2x2_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 2)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 2)





typeA_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_29yr/cond1"
typeB_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_29yr/sim_cond2"

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_90x90_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_90x90_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 90)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 90)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_45x45_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_45x45_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 45)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 45)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_30x30_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_30x30_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 30)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 30)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_18x18_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_18x18_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 18)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 18)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_15x15_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_15x15_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 15)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 15)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_10x10_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_10x10_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 10)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 10)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_9x9_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_9x9_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 9)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 9)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_6x6_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_6x6_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 6)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 6)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_5x5_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_5x5_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 5)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 5)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_3x3_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_3x3_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 3)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 3)

output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_2x2_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_29yr_2x2_filterB_results.csv")

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 2)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 2)
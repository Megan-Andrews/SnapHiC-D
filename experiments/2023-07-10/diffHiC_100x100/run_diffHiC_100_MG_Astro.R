rm(list = ls())
gc() # collect garbage

library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

snap_dir = "/home/maa160/SnapHiC-D/"

source(paste0(snap_dir, "R/utilities.R"))

chrom_size_filepath <- paste0(snap_dir, 'ext/hg19.chrom.sizes')
filter_regions_path <- paste0(snap_dir, "ext/hg19_filter_regions.txt")
TSS_regions_path <- paste0(snap_dir, "/ext/hg19.refGene.transcript.TSS.061421.txt")

experiment_directory <- paste0(snap_dir, "experiments/2023-07-10/diffHiC_100x100/")
result_directory <- paste0(experiment_directory, "results/")

astro_file_list_path <- paste0(snap_dir, "ext/Lee_Astro_samples.txt")
mg_file_list_path <- paste0(snap_dir,"ext/Lee_MG_samples.txt")
typeA_files <- readLines(astro_file_list_path)
typeB_files <- readLines(mg_file_list_path)

chrs = paste0("chr", c(1:22))

print("loading data")
dfs = lapply(c(typeA_files, typeB_files), read_hic_df, chrom_columns=c(2,4),pos_columns=c(3,5), filter_regions_path=filter_regions_path, TSS_regions_path=TSS_regions_path, chrom_size_filepath=chrom_size_filepath, resolution=10000)
# bandnorm_df = bandnorm(hic_df)
# rm(dfs);gc()
# only including valid bins 
for (i in 1:length(dfs)){
  df = dfs[[i]]
  filter_regions = fread(filter_regions_path)
  filter_regions = filter_regions[,c(1:2)]
  colnames(filter_regions) = c('chr','start')
  filter_regions$bin_id = as.integer(filter_regions$start/10000)
  filter_vecs = make_regions_vecs(filter_regions, FALSE, chrom_size_filepath, 
                                  paste0('chr', c(1:22)), 10000)
  is_valid = function(x){filter_vecs[[x['chr']]][as.integer(x['bin1_id'])] & 
      filter_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  valid = apply(df, 1, is_valid)
  valid = as.numeric(valid)
  valid[is.na(valid)] = 0
  
  # # only including TSS regions 
  TSS_regions = fread(TSS_regions_path)
  TSS_regions$bin_id = as.integer((TSS_regions$start+TSS_regions$end)/(2*10000))
  TSS_vecs = make_regions_vecs(TSS_regions, TRUE, chrom_size_filepath, 
                               paste0('chr', c(1:22)), 10000)
  #return(list('TSS_vecs'=TSS_vecs,'df'=binned_df))
  is_TSS = function(x){TSS_vecs[[x['chr']]][as.integer(x['bin1_id'])] |
      TSS_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  TSS = apply(df, 1, is_TSS)
  TSS = as.numeric(TSS)
  TSS[is.na(TSS)] = 0
  dfs[[i]] = df[valid&TSS,]
}

isets = make_iset(dfs)
assayNames(isets) = 'counts'

print("filtering uninteresting bin pairs")  # apply filtering first
keep <- rowSums(assay(isets)>0) > ncol(assay(isets))*0.1  # keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
isets <- isets[keep,]

y <- asDGEList(isets)

# normalization
y <- normOffsets(y, se.out=TRUE)

# modeling and testing 
group = factor(c(rep(1,length(typeA_files)) , rep(2, length(typeB_files))))
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

write.csv(results.r, file = paste0(snap_dir, "experiments/2023-07-10/diffHiC_100x100/results/diffHiC_MG_Astro_results.csv"), row.names = FALSE)


png("Average Abundance Histogram.png")
ave.ab <- aveLogCPM(asDGEList(isets))
hist(ave.ab, xlab="Average abundance", col="grey80", main="")
dev.off()
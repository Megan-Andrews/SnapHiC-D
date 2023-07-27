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

experiment_directory <- paste0(snap_dir, "experiments/2023-07-10/diffHiC_2x2/")
result_directory <- paste0("/project/compbio-lab/scHi-C/Lee2019/results/2023-07-20/")
output_file_filterA <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_30x30_filterA_results.csv")
output_file_filterB <- paste0(result_directory, "diffHiC_Astro_MG_190315_21yr_30x30_filterB_results.csv")

typeA_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_21yr/cond1"
typeB_dir <- "/project/compbio-lab/scHi-C/Lee2019/simulations/Astro_MG_190315_21yr/cond2"

##### BEGIN DIFFHIC #####
print("loading data")

# Read files
run_diffHic <- function(typeA_dir, typeB_dir, output_file, filter_criteria, batch_size){
    dfs = read_hic_cool(typeA_dir, typeB_dir, "chr22")
    isets = make_iset(dfs)
    isets = isets[,c(1:90, (dim(isets)[2]-90+1):dim(isets)[2])] # get first 90 (cond1) and last 90 (cond2) samples
    isets = coarsen_iset(isets, batch_size)

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

run_diffHic(typeA_dir, typeB_dir, output_file_filterA, "A", 30)
run_diffHic(typeA_dir, typeB_dir, output_file_filterB, "B", 30)
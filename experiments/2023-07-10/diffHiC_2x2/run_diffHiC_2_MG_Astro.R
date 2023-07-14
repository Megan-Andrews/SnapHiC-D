rm(list = ls())
gc() # collect garbage

library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

source("/home/maa160/SnapHiC-D/R/utilities.R")

chrom_sizes <- '/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes'

input_directory_Astro <- "/project/compbio-lab/scHi-C/Lee2019/pseudo-bulk_data/Astro/"  
input_directory_MG <- "/project/compbio-lab/scHi-C/Lee2019/pseudo-bulk_data/MG/"  

experiment_directory <- "/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/"
result_directory <- paste0(experiment_directory, "results/")

astro_file_list_path <- paste0(experiment_directory, "file_lists/Astro_10kb_Bulk_file_list.txt") 
mg_file_list_path <- paste0(experiment_directory, "file_lists/MG_10kb_Bulk_file_list.txt") 
astro_file_list <- readLines(astro_file_list_path)
mg_file_list <- readLines(mg_file_list_path)

typeA_files = paste0(input_directory_Astro, astro_file_list)
typeB_files = paste0(input_directory_MG, mg_file_list)

chrs <- paste0("chr", c(1:22, "X", "Y"))

get_diffHiC_results <- function(chr, resolution, typeA_files, typeB_files, chrom_sizes){
    print("loading data")
    diffhic_obj = make_diffhic_object('cool', c(typeA_files,typeB_files), chr,chrom_sizes, resolution)

    print("filtering uninteresting bin pairs") 


    filter_regions <- read.csv("/home/maa160/SnapHiC-D/ext/hg19_filter_regions.txt", sep="\t", header=FALSE)
    gene_transcript <- read.csv("/home/maa160/SnapHiC-D/ext/hg19.refGene.transcript.TSS.061421.txt", sep="\t", header=TRUE)

    colnames(filter_regions) <- c('chr_name', 'x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
    filter_regions <- filter_regions[filter_regions$chr_name == chr, ]
    filter_regions <- filter_regions[, 3] / resolution

    gene_transcript <- gene_transcript[gene_transcript$chr == chr, ]
    gene_transcript$TSS.Bin <- ceiling(gene_transcript$TSS / resolution)
    gene_transcript <- unique(gene_transcript$TSS.Bin)

    # keep <- aveLogCPM(asDGEList(diffhic_obj)) > 2 # should this filtering be applied first
    # diffhic_obj <- diffhic_obj[keep,]

    binIds <- anchorIds(diffhic_obj, type="both")
    keep_filter_regions = !(binIds$first %in% filter_regions) & !(binIds$second %in% filter_regions)
    diffhic_obj <- diffhic_obj[keep_filter_regions,]

    binIds <- anchorIds(diffhic_obj, type="both")
    keep_gene_transcript = (binIds$first %in% gene_transcript) & (binIds$second %in% gene_transcript)
    diffhic_obj <- diffhic_obj[keep_gene_transcript,]
    
    return (diffhic_obj)

}

data = get_diffHiC_results("chr1", 10000, typeA_files, typeB_files, chrom_sizes)     
for (chr in chrs[2:22]){
    print(chr)
    diffhic_obj = get_diffHiC_results(chr, 10000, typeA_files, typeB_files, chrom_sizes)
    print(diffhic_obj)
    data = rbind(data, diffhic_obj)
}
data 

png("Average Abundance Histogram.png")
ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance", col="grey80", main="")
dev.off()

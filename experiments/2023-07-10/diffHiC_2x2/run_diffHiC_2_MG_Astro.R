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
    gene_transcript <- unique(g$TSS.Bin)

    keep <- rowSums(assay(diffhic_obj)>0) > ncol(assay(diffhic_obj))*0.1  # keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
    diffhic_obj <- diffhic_obj[keep,]
    print(keep)
    print(gene_transcript)
    print(filter_regions)

    # i <- interactions(diffhic_obj)
    # rm(keep);gc()

    # print("normalization")
    # y <- normOffsets(asDGEList(diffhic_obj), se.out=TRUE)

    # rm(diffhic_obj);gc()
    # return([y, i])    
}

get_diffHiC_results("chr21", 10000, typeA_files, typeB_files, chrom_sizes)
# print("modeling and testing" )
# group = factor(c(rep(1, length(typeA_files)*length(chrs)), rep(2, length(typeB_files)*length(chrs))))
# design = model.matrix(~group)
# rm(group);gc()

# all_DGE_List <- NULL
# interactions <- NULL
# for chr in chrs{
#     result = get_diffHiC_results(chr)
#     y = result[1]
#     i = result[2]
#     if (is.null(all_DGE_List)){
#         all_DGE_List = y
#     }else{  
#         all_DGE_List = rbind(all_DGE_List, y)
#     }
#     if (is.null(interactions)){
#         interactions = i
#     }else{
#         interactions = rbind(interactions, i)
#     }
# }


# y <- estimateDisp(y, design) # From edgeR 
# print(object_size(y)) # 809.80 MB

# gc()
# print("glmQLFit")
# fit <- glmQLFit(y, design, robust=TRUE)
# rm(design)
# rm(y)
# gc()

# print("glmQLFTest")
# result <- glmQLFTest(fit)
# print(object_size(result))

# gc()
# print("adj.p")
# adj.p <- p.adjust(result$table$PValue, method="BH")
# sum(adj.p <= 0.05)
# print(object_size(adj.p))

# print("useful.cols")
# useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
# inter.frame <- as.data.frame(interactions)[,useful.cols]
# results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
# o.r <- order(results.r$PValue)
# results.r = results.r[o.r,]
# rm(results.r, useful.cols, inter.frame, adj.p, o.r);gc()


# write.csv(all_chr_result, file = paste0(result_directory, "diffHiC_MG_Astro_results.csv") , row.names = FALSE)













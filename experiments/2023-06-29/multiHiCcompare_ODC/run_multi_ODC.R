library(edgeR)
library(diffHic)
library(csaw)
library(statmod)
library(dplyr)

source("/home/maa160/SnapHiC-D/R/utilities.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"  
A_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-29/multiHiCcompare_ODC/file_lists/ODC_1_100kb_file_list.txt"  
B_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-29/multiHiCcompare_ODC/file_lists/ODC_2_100kb_file_list.txt" 
A_file_list <- readLines(A_file_list_path)
B_file_list <- readLines(B_file_list_path)

A_num = length(A_file_list) 
B_num = length(B_file_list) 

typeA_files = paste0(input_directory, A_file_list)
typeB_files = paste0(input_directory, B_file_list)


chrom_sizes = '/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes'


### multiHiCcompare 

# loading data (intrinsic filtering)
hicexp = make_multiHiCcompare_object('txt', c(typeA_files,typeB_files), 'chr22',
                                     chrom_sizes, 100000,
                                     c(rep('A', 449), rep('B', 422)))
# normalization 
MD_hicexp(hicexp,plot.loess = TRUE)
norm_hicexp = cyclic_loess(hicexp, verbose = FALSE, 
                           parallel = FALSE, span = 0.2)
MD_hicexp(norm_hicexp,plot.loess = TRUE)
# modeling and testing 
norm_hicexp = hic_exactTest(norm_hicexp, p.method = 'fdr')
norm_hicexp = hic_glm(norm_hicexp,design,coef=ncol(design))


# analysis 
diffhic.res = results.r[,c('start1', 'start2', 'logFC', 'PValue', 'FDR')]
multi.res = results(norm_hicexp)[,c('region1', 'region2', 'logFC', 'p.value', 'p.adj')]
colnames(diffhic.res) = c('region1', 'region2', 'logFC', 'p.value', 'p.adj')
diffhic.res[,c('region2','region1')] = diffhic.res[,c('region1','region2')]*100000
res = inner_join(x=diffhic.res,multi.res,by=c('region1'='region1', 'region2'='region2'))


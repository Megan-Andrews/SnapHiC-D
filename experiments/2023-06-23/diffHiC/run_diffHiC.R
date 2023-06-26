library(diffHic)
library(rhdf5)
library(GenomicRanges)
library(InteractionSet)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)
source("/home/maa160/SnapHiC-D/R/utilities.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"  #should we run diffHiC on 10kb data 100kb data or the imputed data? 
astro_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/Astro_100kb_file_list.txt"  
mg_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/MG_100kb_file_list.txt" 
astro_file_list <- readLines(astro_file_list_path)
mg_file_list <- readLines(mg_file_list_path)

typeA_files = paste0(input_directory, astro_file_list[1:30])
typeB_files = paste0(input_directory, mg_file_list[1:30])

# loading data
diffhic_obj = make_diffhic_object('cool', c(typeA_files,typeB_files), 'chr22','/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes', 100000)





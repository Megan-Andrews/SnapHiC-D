rm(list = ls())
gc() # collect garbage

library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

snap_dir = "/Users/megan/Projects/SnapHiC-D/"

source(paste0(snap_dir, "R/utilities.R"))

chrom_size_filepath <- paste0(snap_dir, 'ext/hg19.chrom.sizes')
filter_regions_path <- paste0(snap_dir, "ext/hg19_filter_regions.txt")
TSS_regions_path <- paste0(snap_dir, "/ext/hg19.refGene.transcript.TSS.061421.txt")

experiment_directory <- paste0(snap_dir, "experiments/2023-07-10/diffHiC_2x2/")
result_directory <- paste0(experiment_directory, "results/")
output_directory <- paste0(result_directory, "diffHiC_MG_Astro_results.csv")

#astro_file_list_path <- paste0(snap_dir, "ext/Lee_Astro_samples.txt")
#mg_file_list_path <- paste0(snap_dir,"ext/Lee_MG_samples.txt")
#typeA_files <- readLines(astro_file_list_path)
#typeB_files <- readLines(mg_file_list_path)

typeA_files <- c(
  "GSM3748102_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD001_indexed_contacts.txt.gz",
  "GSM3748103_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD004_indexed_contacts.txt.gz",
  "GSM3748104_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD006_indexed_contacts.txt.gz",
  "GSM3748105_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD007_indexed_contacts.txt.gz",
  "GSM3748106_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD008_indexed_contacts.txt.gz",
  "GSM3748102_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD001_indexed_contacts.txt.gz",
  "GSM3748103_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD004_indexed_contacts.txt.gz",
  "GSM3748104_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD006_indexed_contacts.txt.gz",
  "GSM3748105_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD007_indexed_contacts.txt.gz",
  "GSM3748106_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A10_AD008_indexed_contacts.txt.gz"
)
typeA_files <- paste0("/Users/megan/CompBio-Data/GSE130711_RAW/", typeA_files)
typeB_files <- c(
  "GSM3748121_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD012_indexed_contacts.txt.gz",
  "GSM3748120_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD010_indexed_contacts.txt.gz",
  "GSM3748119_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD008_indexed_contacts.txt.gz",
  "GSM3748118_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD007_indexed_contacts.txt.gz",
  "GSM3748117_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD006_indexed_contacts.txt.gz",
  "GSM3748121_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD012_indexed_contacts.txt.gz",
  "GSM3748120_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD010_indexed_contacts.txt.gz",
  "GSM3748119_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD008_indexed_contacts.txt.gz",
  "GSM3748118_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD007_indexed_contacts.txt.gz",
  "GSM3748117_181218_snm3Cseq_hs_21yr_BA10_UMB5577_1_UMB5577_2_A12_AD006_indexed_contacts.txt.gz"
)
typeB_files <- paste0("/Users/megan/CompBio-Data/GSE130711_RAW/", typeB_files)

applyBandNorm = FALSE




##### BEGIN DIFFHIC #####
print("loading data")

# Read typeA files
dfs_A = lapply(c(typeA_files), read_hic_df, chrom_columns=c(2,4),pos_columns=c(3,5), filter_regions_path=filter_regions_path, TSS_regions_path=TSS_regions_path, chrom_size_filepath=chrom_size_filepath, resolution=10000)

# Read typeB files
dfs_B = lapply(c(typeB_files), read_hic_df, chrom_columns=c(2,4),pos_columns=c(3,5), filter_regions_path=filter_regions_path, TSS_regions_path=TSS_regions_path, chrom_size_filepath=chrom_size_filepath, resolution=10000)

if (applyBandNorm==TRUE){
  print("applying BandNorm")
  dfs_A = bandnorm(rbind_hic_dfs(dfs_A))
  dfs_A = lapply(array(split(dfs_A, dfs_A$cell)), function(x) x[,c("chr", "bin1_id", "bin2_id", "BandNorm")])
  
  dfs_B = bandnorm(rbind_hic_dfs(dfs_B))
  dfs_B = lapply(array(split(dfs_B, dfs_B$cell)), function(x) x[,c("chr", "bin1_id", "bin2_id", "BandNorm")])
}

dfs = c(dfs_A, dfs_B)
rm(dfs_A,dfs_B);gc()

# only including valid bins 
# dfs = lapply(dfs, filter_df, filter_regions_path=filter_regions_path, TSS_regions_path=TSS_regions_path)

isets = make_iset(dfs)

isets = coarsen_iset(isets, 5)

assayNames(isets) = 'counts'
anchors = anchors(isets,type="both")
GI = GInteractions(anchors$first, anchors$second)

y <- asDGEList(isets)
y <- normOffsets(y, se.out=TRUE) # normalization
counts = y$counts
rownames(counts) <- NULL
colnames(counts) <- NULL

chrs = elementMetadata(isets)$X

isets = InteractionSet(as(counts, "matrix"), GI)
elementMetadata(isets) = chrs

isets = filter_iset(isets, filter_regions_path, TSS_regions_path)

print("filtering uninteresting bin pairs")  # apply filtering first
keep <- rowSums(assay(isets)>0) > ncol(assay(isets))*0.1  # keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
isets <- isets[keep,]


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

write.csv(results.r, file = paste0(snap_dir, output_directory), row.names = FALSE)

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

    keep <- aveLogCPM(asDGEList(diffhic_obj)) > 2 # should this filtering be applied first
    diffhic_obj <- diffhic_obj[keep,]

    binIds <- anchorIds(diffhic_obj, type="both")
    keep_filter_regions = !(binIds$first %in% filter_regions) & !(binIds$second %in% filter_regions)
    diffhic_obj <- diffhic_obj[keep_filter_regions,]

    binIds <- anchorIds(diffhic_obj, type="both")
    keep_gene_transcript = (binIds$first %in% gene_transcript) & (binIds$second %in% gene_transcript)
    diffhic_obj <- diffhic_obj[keep_gene_transcript,]
    
    return (diffhic_obj)

}

data = get_diffHiC_results("chr17", 10000, typeA_files, typeB_files, chrom_sizes)     
for (chr in chrs[18:22]){
    print(chr)
    diffhic_obj = get_diffHiC_results(chr, 10000, typeA_files, typeB_files, chrom_sizes)
    print(diffhic_obj)
    data = rbind(data, diffhic_obj)
    rm(diffhic_obj); gc()
}
data 

y <- asDGEList(data)

# normalization
y <- normOffsets(y, se.out=TRUE)

# par(mfrow=c(1,2))
# ab <- aveLogCPM(asDGEList(data))
# o <- order(ab)
# adj.counts <- cpm(asDGEList(data), log=TRUE)
# mval <- adj.counts[,3]-adj.counts[,2]

# smoothScatter(ab, mval, xlab="A", ylab="M", 
#               main="MG (1) vs. Astro (2) \n before normalization")
# fit <- loessFit(x=ab, y=mval)
# lines(ab[o], fit$fitted[o], col="red")

# savePlot(filename = "before_normalization.png", type = "png", width = 8, height = 6)


# ab <- aveLogCPM(y)
# o <- order(ab)
# adj.counts <- cpm(y, log=TRUE)
# mval <- adj.counts[,3]-adj.counts[,2]
# smoothScatter(ab, mval, xlab="A", ylab="M", 
#               main="MG (1) vs. Astro (2) \n after normalization")
# fit <- loessFit(x=ab, y=mval)
# lines(ab[o], fit$fitted[o], col="red")
# print("saving")
# savePlot(filename = "after_normalization.png", type = "png", width = 8, height = 6)

# dev.off()


# modeling and testing 
group = factor(c(rep(1,length(typeA_files)) , rep(2, length(typeB_files))))
design = model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
result <- glmQLFTest(fit)
adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= 0.05)
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
inter.frame <- as.data.frame(interactions(data))[,useful.cols]
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
o.r <- order(results.r$PValue)
results.r = results.r[o.r,]

write.csv(results.r, file = "/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/diffHiC_MG_Astro_results.csv", row.names = FALSE)


# png("Average Abundance Histogram.png")
# ave.ab <- aveLogCPM(asDGEList(data))
# hist(ave.ab, xlab="Average abundance", col="grey80", main="")
# dev.off()

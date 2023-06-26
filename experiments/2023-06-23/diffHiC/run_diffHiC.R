library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

source("/home/maa160/SnapHiC-D/R/utilities.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"  
astro_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/Astro_100kb_file_list.txt"  
mg_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/MG_100kb_file_list.txt" 
astro_file_list <- readLines(astro_file_list_path)
mg_file_list <- readLines(mg_file_list_path)

typeA_files = paste0(input_directory, astro_file_list)
typeB_files = paste0(input_directory, mg_file_list)

print("loading data")
diffhic_obj = make_diffhic_object('cool', c(typeA_files,typeB_files), 'chr22','/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes', 100000)

print("filtering uninteresting bin pairs") 
keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
diffhic_obj <- diffhic_obj[keep,]
y <- asDGEList(diffhic_obj)

print("normalization")
y <- normOffsets(y, se.out=TRUE)

print("visualization")
par(mfrow=c(1,2))
png(filename = "loess_smoothing.png", width = 800, height = 600)
ab <- aveLogCPM(asDGEList(diffhic_obj))
o <- order(ab)

adj.counts <- cpm(asDGEList(diffhic_obj), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="Astro (1) vs. MG (2) \n before normalization")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")


ab <- aveLogCPM(y)
o <- order(ab)
adj.counts <- cpm(y, log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="Astro (1) vs. MG (2) \n after normalization")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
dev.off()

print("modeling and testing" )
group = factor(c(rep(1, 449), rep(2, 422)))
design = model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
result <- glmQLFTest(fit)
adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= 0.05)
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
inter.frame <- as.data.frame(interactions(diffhic_obj))[,useful.cols]
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
o.r <- order(results.r$PValue)
results.r = results.r[o.r,]
write.csv(results.r, file = "diffHiC_results.csv", row.names = FALSE)





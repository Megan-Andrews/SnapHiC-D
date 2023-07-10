rm(list = ls())
gc() # collect garbage

library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

source("/home/maa160/SnapHiC-D/R/utilities.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/"  
result_directory <- "/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC/results/"
astro_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC/file_lists/Astro_10kb_file_list.txt"  
mg_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC/file_lists/MG_10kb_file_list.txt" 
astro_file_list <- readLines(astro_file_list_path)
mg_file_list <- readLines(mg_file_list_path)

typeA_files = paste0(input_directory, astro_file_list)
typeB_files = paste0(input_directory, mg_file_list)

print("loading data")
diffhic_obj = make_diffhic_object('cool', c(typeA_files,typeB_files), 'chr22','/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes', 100000)
rm(cool2matrix, cool2sparse, df2mat, make_diffhic_object, make_multiHiCcompare_object)
rm(input_directory, astro_file_list_path, mg_file_list_path, astro_file_list, mg_file_list,typeB_files, typeA_files)
gc()

print("filtering uninteresting bin pairs") 
keep <- rowSums(assay(diffhic_obj)>0) > ncol(assay(diffhic_obj))*0.1  # keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
diffhic_obj <- diffhic_obj[keep,]
y <- asDGEList(diffhic_obj)
i <- interactions(diffhic_obj)
print(object_size(y))
rm(keep)
gc()

print("normalization")
y <- normOffsets(asDGEList(diffhic_obj), se.out=TRUE)
print(paste("Object size: y ", object_size(y)))


print("visualization")
# par(mfrow=c(1,2))
png(filename = paste0(result_directory, "before_loess_smoothing.png"), width = 800, height = 600)


print("aveLogCPM")
print(object_size(y))
ab <- aveLogCPM(asDGEList(diffhic_obj))
print(object_size(ab))
o <- order(ab)

print("adjusted counts")
adj.counts <- cpm(asDGEList(diffhic_obj) , log=TRUE)
rm(diffhic_obj)
gc()

print(object_size(adj.counts ))
mval <- adj.counts[,3]-adj.counts[,2]

print("smoothScatter")
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="Astro (1) vs. MG (2) \n before normalization")

print("loessFit")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

dev.off()
png(filename = paste0(result_directory, "after_loess_smoothing.png") , width = 800, height = 600)
gc()

print("aveLogCPM")
ab <- aveLogCPM(y)
print(object_size(ab))
o <- order(ab)
gc()

print("adj.counts")
adj.counts <- cpm(y, log=TRUE)
print(object_size(adj.counts )) # 404.16 MB
mval <- adj.counts[,3]-adj.counts[,2]
gc()

print("smoothScatter")
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="Astro (1) vs. MG (2) \n after normalization")

fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
dev.off()

rm(fit)
rm(ab)
rm(o)
rm(mval)
rm(adj.counts)
gc()

print("modeling and testing" )
group = factor(c(rep(1, 449), rep(2, 422)))
design = model.matrix(~group)
rm(group)

y <- estimateDisp(y, design) # From edgeR 
print(object_size(y)) # 809.80 MB

gc()
print("glmQLFit")
fit <- glmQLFit(y, design, robust=TRUE)
rm(design)
rm(y)
gc()

print("glmQLFTest")
result <- glmQLFTest(fit)
print(object_size(result))

gc()
print("adj.p")
adj.p <- p.adjust(result$table$PValue, method="BH")
sum(adj.p <= 0.05)
print(object_size(adj.p))

print("useful.cols")
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
print(object_size(adj.p))
inter.frame <- as.data.frame(i)[,useful.cols]
print(object_size(inter.frame))
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
print(object_size(results.r ))
o.r <- order(results.r$PValue)
results.r = results.r[o.r,]
print(object_size(results.r ))

write.csv(results.r, file = paste0(result_directory, "diffHiC_MG_Astro_results.csv") , row.names = FALSE)





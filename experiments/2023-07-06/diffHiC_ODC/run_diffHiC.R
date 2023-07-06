rm(list = ls())
gc()
library(diffHic)
library(rhdf5)
library(edgeR)
library(csaw)
library(statmod)
library(dplyr)

source("/home/maa160/SnapHiC-D/R/utilities.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"  
A_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-07-06/diffHiC_ODC/file_lists/ODC_1_100kb_file_list.txt"  
B_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-07-06/diffHiC_ODC/file_lists/ODC_2_100kb_file_list.txt" 
A_file_list <- readLines(A_file_list_path)
B_file_list <- readLines(B_file_list_path)

typeA_files = paste0(input_directory, A_file_list)
typeB_files = paste0(input_directory, B_file_list)

A_num = length(A_file_list)
B_num = length(B_file_list)

chrom_sizes = '/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes'
output_path = "/home/maa160/SnapHiC-D/experiments/2023-07-06/diffHiC_ODC/results/"

print("loading data")
diffhic_obj = make_diffhic_object('cool', c(typeA_files,typeB_files), 'chr22',chrom_sizes, 100000)
rm(cool2matrix, cool2sparse, df2mat, make_diffhic_object, make_multiHiCcompare_object)
rm(input_directory, astro_file_list_path, mg_file_list_path, astro_file_list)
rm(mg_file_list,typeB_files, typeA_files)

gc()
print("filtering uninteresting bin pairs") 
# keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
keep <- rowSums(assay(diffhic_obj)>0) > ncol(assay(diffhic_obj))*0.1
diffhic_obj <- diffhic_obj[keep,]
y <- asDGEList(diffhic_obj)
i <- interactions(diffhic_obj)
print(object_size(y))
rm(keep)

gc()
print("normalization")
y <- normOffsets(asDGEList(diffhic_obj), se.out=TRUE)
print(object_size(y))

gc()
print("visualization")
# par(mfrow=c(1,2))
png(filename = paste0(output_path, "before_loess_smoothing.png"), width = 800, height = 600)

gc()
print("aveLogCPM")
print(object_size(y))
ab <- aveLogCPM(asDGEList(diffhic_obj))
print(object_size(ab))
o <- order(ab)

gc()
print("adjusted counts")
adj.counts <- cpm(asDGEList(diffhic_obj) , log=TRUE)
rm(diffhic_obj)
gc()

print(object_size(adj.counts ))
mval <- adj.counts[,3]-adj.counts[,2]

gc()
print("smoothScatter")
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="ODC (1) vs. ODC (2) \n before normalization")

gc()
print("loessFit")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

dev.off()
png(filename = paste0(output_path, "after_loess_smoothing.png"), width = 800, height = 600)

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
              main="ODC (1) vs. ODC (2) \n after normalization")

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
group = factor(c(rep(1, A_num), rep(2, B_num)))
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

write.csv(results.r, file = paste0(output_path, "diffHiC_MG_Astro_results.csv"), row.names = FALSE)





source('R/utilities.R')
library(edgeR)
library(diffHic)
library(csaw)
library(statmod)
library(dplyr)


typeA_files = c('/media/nedooshki/drive1/VSS-Hi-C/data/GM12878/Hi-C/rep1/observed_raw_100Kb/chr21_chr21.txt',
                '/media/nedooshki/drive1/VSS-Hi-C/data/GM12878/Hi-C/rep2/observed_raw_100Kb/chr21_chr21.txt')

typeB_files = c('/media/nedooshki/drive1/VSS-Hi-C/data/K562/Hi-C/rep1/observed_raw_100Kb/chr21_chr21.txt',
                '/media/nedooshki/drive1/VSS-Hi-C/data/K562/Hi-C/rep2/observed_raw_100Kb/chr21_chr21.txt')

### diffHiC

# loading data
diffhic_obj = make_diffhic_object('txt', c(typeA_files,typeB_files), 'chr21',
                           '/media/nedooshki/drive1/IChDA/data/supp/hg19.chrom.sizes', 100000)
# filtering uninteresting bin pairs 
keep <- aveLogCPM(asDGEList(diffhic_obj)) > 0
diffhic_obj <- diffhic_obj[keep,]
y <- asDGEList(diffhic_obj)
# normalization
y <- normOffsets(y, se.out=TRUE)

par(mfrow=c(1,2))
ab <- aveLogCPM(asDGEList(diffhic_obj))
o <- order(ab)
adj.counts <- cpm(asDGEList(diffhic_obj), log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="K562 (1) vs. GM12878 (2) \n before normalization")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

ab <- aveLogCPM(y)
o <- order(ab)
adj.counts <- cpm(y, log=TRUE)
mval <- adj.counts[,3]-adj.counts[,2]
smoothScatter(ab, mval, xlab="A", ylab="M", 
              main="K562 (1) vs. GM12878 (2) \n after normalization")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")
dev.off()
# modeling and testing 
group = factor(c(1,1,2,2))
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



### multiHiCcompare 

# loading data (intrinsic filtering)
hicexp = make_multiHiCcompare_object('txt', c(typeA_files,typeB_files), 'chr21',
                                     '/media/nedooshki/drive1/IChDA/data/supp/hg19.chrom.sizes', 100000,
                                     c('A','A','B','B'))
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


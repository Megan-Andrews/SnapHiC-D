library(ggplot2)
library(ggpubr)
library(edgeR)
source("R/utilities.R")


# load data 
Astro_dir = "/media/nedooshki/drive1/scHi-C data/Lee/Astro/"
MG_dir = "/media/nedooshki/drive1/scHi-C data/Lee/MG/"
Astro_files = list.files(Astro_dir)
MG_files = list.files(MG_dir)
Astro_filepaths = paste0(Astro_dir, Astro_files)
MG_filepaths = paste0(MG_dir, MG_files)

Astro_data = lapply(Astro_filepaths, read_hic_df, c(2,4),c(3,5),100000,
                    'ext/hg19_filter_regions.txt', 
                    'ext/hg19.refGene.transcript.TSS.061421.txt', 'ext/hg19.chrom.sizes')

MG_data = lapply(MG_filepaths, read_hic_df, c(2,4),c(3,5),100000,
                    'ext/hg19_filter_regions.txt', 
                    'ext/hg19.refGene.transcript.TSS.061421.txt', 'ext/hg19.chrom.sizes')

all_df = rbind_hic_dfs(c(Astro_data,MG_data))

groups = c(rep('Astro',100),rep('MG',100))

# raw 
raw_embedding = hic_embedding(all_df, 'df', groups)
raw_plot = ggplot(raw_embedding, aes(x = X1, y = X2, col=cell_type)) + geom_point() + ylab("UMAP 2") +
  xlab("UMAP 1") + theme_bw(base_size = 15) + ggtitle('raw')



# BandNorm 
bandnorm_df <- bandnorm(all_df)
colnames(bandnorm_df) = c('chr', 'bin1_id', 'bin2_id', 'cell', 'diag', 'count')
bandnorm_embedding = hic_embedding(bandnorm_df, 'df', groups)
bandnorm_plot = ggplot(bandnorm_embedding, aes(x = X1, y = X2, col=cell_type)) + geom_point() + ylab("UMAP 2") +
  xlab("UMAP 1") + theme_bw(base_size = 15) + ggtitle('BandNorm')



# multiHiCcompare (fastlo on pools)
hicexp = make_hicexp_obj(c(Astro_data,MG_data), 100000, groups)
fastlo_hicexp = fastlo(hicexp)
fastlo_mat = t(data.frame(hic_table(fastlo_hicexp))[,5:204])
fastlo_embedding = hic_embedding(fastlo_mat, 'mat', groups)
fastlo_plot = ggplot(fastlo_embedding, aes(x = X1, y = X2, col=cell_type)) + geom_point() + ylab("UMAP 2") +
  xlab("UMAP 1") + theme_bw(base_size = 15) + ggtitle('fastlo on pools')


# diffHiC (fastlo on MA plot)
IFs = hic_table(hicexp)[,5:204]
edgeR_obj = DGEList(IFs)
edgeR_obj = normOffsets(edgeR_obj, se.out=TRUE)
cpm_IFs = cpm(edgeR_obj, log=FALSE)
cpm_embedding = hic_embedding(t(cpm_IFs), 'mat', groups)
cpm_plot = ggplot(cpm_embedding, aes(x = X1, y = X2, col=cell_type)) + geom_point() + ylab("UMAP 2") +
  xlab("UMAP 1") + theme_bw(base_size = 15) + ggtitle('fastlo on MA')


plt <- ggarrange(raw_plot, bandnorm_plot, fastlo_plot, cpm_plot,
          ncol = 2, nrow = 2, common.legend = TRUE)
ggsave("results/2023-07-12/normalization&clustering-wg.png", plt,
       width = 20, height = 20, units = "cm", bg = "white")

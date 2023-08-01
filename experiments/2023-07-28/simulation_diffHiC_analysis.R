source("R/utilities.R")
source("R/evaluation_utilities.R")

filter_regions_path = 'ext/hg19_filter_regions.txt'
ground_truth_DCC_path = "/media/nedooshki/drive1/scHi-C data/Lee2019/simulation/Astro_MG_190315_21yr_100kb/DCC.txt"
typeA_dir <- "/media/nedooshki/drive1/scHi-C data/Lee2019/simulation/Astro_MG_190315_21yr_100kb/sim_cond1/"
typeB_dir <- "/media/nedooshki/drive1/scHi-C data/Lee2019/simulation/Astro_MG_190315_21yr_100kb/sim_cond2/"

dfs = read_hic_cool(typeA_dir, typeB_dir, "chr22")
DCCs = read.table(ground_truth_DCC_path, col.names = c('bin1_id', 'bin2_id', 'LFC'))

experiment_sizes = c(2, 10, 15, 30, 45, 90)
isets = prep_isets(dfs, filter_regions_path, experiment_sizes)

results = list()
for (es in c(2, 10, 15, 30, 45, 90)){
  results[[es]] = make_diffHiC_results(isets[[es]], es)
}



merged_results = list()
for (es in c(2, 10, 15, 30, 45, 90)){
  merged_results[[es]] = merge(results[[es]]$result, DCCs, by.x = c('start1', 'start2'), 
            by.y = c('bin1_id', 'bin2_id'), all = TRUE)
  merged_results[[es]]$LFC[is.na(merged_results[[es]]$LFC)] = 0
  merged_results[[es]]$logFC[is.na(merged_results[[es]]$logFC)] = 0
  merged_results[[es]]$FDR[is.na(merged_results[[es]]$FDR)] = 1
}


# pvalue vs LFC plot

volcano_plots <- lapply(experiment_sizes, function(es){
  ggplot(results[[es]]$result) + 
    geom_point(aes(logFC, -log10(PValue))) + 
    labs(x = "logFC", y = "-log10(p-value)", title = paste0("experiment size: ", es))
})

volcano_plt = ggarrange(volcano_plots[[1]], volcano_plots[[2]], volcano_plots[[3]],
                        volcano_plots[[4]], volcano_plots[[5]], volcano_plots[[6]], 
                        nrow = 2, ncol = 3)
ggsave("results/2023-07-28/simulation_LFC_pval_volcano.png", volcano_plt,
       width = 20, height = 11, units = "cm", bg = "white")

# confusion scatter plot

plots = lapply(experiment_sizes, function(es){
  get_confusion_scatter_plt(merged_results[[es]], es, 0.05)
})
plt = ggarrange(plots[[1]], plots[[2]], plots[[3]],
                plots[[4]], plots[[5]], plots[[6]], 
                nrow = 2, ncol = 3, common.legend = TRUE) 
ggsave("results/2023-07-28/pseudobulking_results.png", plt,
       width = 20, height = 11, units = "cm", bg = "white")

# Accuracy, precision, recall 

stats_df = NULL
for (cfdr in c(0.01, 0.025, 0.05, 0.075, 0.1)){
  for (es in c(2, 10, 15, 30, 45, 90)){
    stats = calculate_stats(merged_results[[es]], cfdr)
    if (is.null(acc_df)){
      stats_df= data.frame(es = es, cfdr = cfdr, precision = stats[['precision']],
                           recall = stats[['recall']], accuracy = stats[['accuracy']])
    }
    else{
      stats_df = rbind(stats_df, data.frame(es = es, cfdr = cfdr, precision = stats[['precision']],
                                            recall = stats[['recall']], 
                                            accuracy = stats[['accuracy']]))
    }
  }
}
stats_df$es = as.factor(stats_df$es)

accuracy_plt = ggplot(data = stats_df) + 
  geom_line(aes(cfdr, accuracy, colour = es)) +
  labs(x = 'Target FDR', y = 'Accuracy') +
  scale_color_discrete(name = "Experiment size")

recall_plt = ggplot(data = stats_df) + 
  geom_line(aes(cfdr, recall, colour = es)) +
  labs(x = 'Target FDR', y = 'Recall') +
  scale_color_discrete(name = "Experiment size")

precision_plt = ggplot(data = stats_df) + 
  geom_line(aes(cfdr, precision, colour = es)) +
  labs(x = 'Target FDR', y = 'Precision') +
  scale_color_discrete(name = "Experiment size")

ggsave("results/2023-07-28/pseudobulking&accuracy.png", accuracy_plt,
       width = 17, height = 10, units = "cm", bg = "white")
ggsave("results/2023-07-28/pseudobulking&precision.png", precision_plt,
       width = 17, height = 10, units = "cm", bg = "white")
ggsave("results/2023-07-28/pseudobulking&recall.png", recall_plt,
       width = 17, height = 10, units = "cm", bg = "white")

# controlling FDR 

fdr_df = NULL
for (cfdr in c(0.01, 0.025, 0.05, 0.075, 0.1)){
  for (es in c(2, 10, 15, 30, 45, 90)){
    observed_fdr = get_observed_FDR(merged_results[[es]], cfdr)
    if (is.null(fdr_df)){
      fdr_df= data.frame(es = es, cfdr = cfdr, observed_fdr = observed_fdr)
    }
    else{
      fdr_df = rbind(fdr_df, data.frame(es = es, cfdr = cfdr, observed_fdr = observed_fdr))
    }
  }
}
fdr_df$es = as.factor(fdr_df$es)

fdr_plt = ggplot(data = fdr_df) + 
  geom_line(aes(cfdr, observed_fdr, colour = es)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  labs(x = 'Target FDR', y = 'Observed FDR') +
  scale_color_discrete(name = "Experiment size")

ggsave("results/2023-07-28/pseudobulking&FDR.png", fdr_plt,
       width = 17, height = 10, units = "cm", bg = "white")

# LFC error 

LFC_errors = NULL
for (es in experiment_sizes){
  if (is.null(LFC_errors)){
    LFC_errors = data.frame(es = es, 
                            LFC_error = (merged_results[[es]]$logFC - 
                                           merged_results[[es]]$LFC)^2)
  }
  else{
    LFC_errors = rbind(LFC_errors, 
                       data.frame(es = es, 
                                  LFC_error = (merged_results[[es]]$logFC - 
                                                 merged_results[[es]]$LFC)^2))
  }
}
LFC_errors$es = as.factor(LFC_errors$es)
LFC_error_plt = ggplot(LFC_errors, aes(es, LFC_error, colour = es)) + 
  geom_boxplot() + labs(x = "Experiment size", y = "LFC error") + 
  theme(legend.position = "none")
ylim1 = boxplot.stats(LFC_errors$LFC_error)$stats[c(1, 5)]
LFC_error_plt = LFC_error_plt + coord_cartesian(ylim = ylim1*1.05)

ggsave("results/2023-07-28/pseudobulking&LFC_error.png", LFC_error_plt,
       width = 13, height = 8, units = "cm", bg = "white")

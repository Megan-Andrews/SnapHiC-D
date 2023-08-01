library(ggplot2)

get_confusion_scatter_plt <- function(merged_res, es, control_FDR){
  is_FDR_DCC = merged_res$FDR < control_FDR 
  is_gt_DCC = abs(merged_res$LFC) > 0 
  bp_type = vector(mode = "character", length = length(is_FDR_DCC))
  bp_type[is_FDR_DCC & is_gt_DCC] = 'TP'
  bp_type[!is_FDR_DCC & !is_gt_DCC] = 'TN'
  bp_type[is_FDR_DCC & !is_gt_DCC] = 'FP'
  bp_type[!is_FDR_DCC & is_gt_DCC] = 'FN'
  merged_res$type = bp_type
  experiment_name = paste0("diffHiC (", es, "*", es, ")")
  plt = ggplot(merged_res) + 
    geom_point(aes(x = LFC, y = logFC, colour = type)) +
    labs(x = "ground truth LFC", y = paste0(experiment_name, " LFC")) + 
    scale_color_discrete(name = "")
  plt
}

calculate_stats <- function(merged_res, control_FDR){
  is_FDR_DCC = merged_res$FDR < control_FDR 
  is_gt_DCC = abs(merged_res$LFC) > 0 
  TP = sum(is_FDR_DCC & is_gt_DCC)
  TN = sum(!is_FDR_DCC & !is_gt_DCC)
  FP = sum(is_FDR_DCC & !is_gt_DCC)
  FN = sum(!is_FDR_DCC & is_gt_DCC)
  #print(paste0("TP: ", TP, ", TN: ", TN, ", FP: ", FP, ", FN: ", FN))
  list('recall' = TP/(TP+FN), 'precision' = TP/(TP+FP), 
       'accuracy' = (TP + TN) / (TP + TN + FP + FN))
}

get_observed_FDR <- function(merged_res, control_FDR){
  FDR_DCCs = merged_res[merged_res$FDR < control_FDR,]
  null_DCCs = FDR_DCCs[FDR_DCCs$LFC==0,]
  #print(dim(FDR_DCCs))
  #print(dim(null_DCCs))
  dim(null_DCCs)[1]/dim(FDR_DCCs)[1]
}
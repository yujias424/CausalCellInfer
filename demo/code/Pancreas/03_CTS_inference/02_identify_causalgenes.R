source("/mnt/md0/yujia/project/github_package/CausalCellInfer/bin/cell_marker_identification/secondary_analyses.R")
library(data.table)

remove0sd <- function(dataframe){
  dataframe[,which(apply(dataframe, 2, function(x) length(unique(x)))>1)]
}

tmp_cell_exp <- as.data.frame(data.table::fread("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/res/Pancreas/CTS/alpha_cell.csv.gz"))
rownames(tmp_cell_exp) <- tmp_cell_exp$V1
tmp_cell_exp$V1 <- NULL
tmp_cell_exp <- t(tmp_cell_exp)
tmp_cell_exp <- remove0sd(tmp_cell_exp)

tmp_cell_exp[1:10, 1:10]

meta <- as.data.frame(fread("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/bulk/Pancreas/GSE50244_meta.csv"))
meta <- meta[!is.na(meta$Hba1c), ]
meta_control <- meta[meta$Hba1c < 6, ]
meta_control$T2D <- 0
meta_case <- meta[meta$Hba1c >= 6.5, ]
meta_case$T2D <- 1
meta <- rbind(meta_control, meta_case)
meta$ID <- paste0("Sample_", meta$ID)

causal_genes <- get_causal_genes(tmp_cell_exp, meta$T2D, 0.05, 10)

causal_genes


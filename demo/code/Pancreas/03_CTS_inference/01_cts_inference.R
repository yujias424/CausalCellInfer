library(ENIGMA)
library(data.table)
library(zellkonverter)
library(scater)
library(abind)

# Load the Bulk data
bulk <- as.data.frame(fread("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz"))
colnames(bulk) <- bulk[1, ]
bulk <- bulk[-1, ]
row.names(bulk) <- bulk$id
bulk$id <- NULL
colnames(bulk) <- paste0("Sample_", colnames(bulk))

meta <- as.data.frame(fread("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_meta.csv"))
meta <- meta[!is.na(meta$Hba1c), ]
meta_control <- meta[meta$Hba1c < 6, ]
meta_control$T2D <- 0
meta_case <- meta[meta$Hba1c >= 6.5, ]
meta_case$T2D <- 1
meta <- rbind(meta_control, meta_case)
meta$ID <- paste0("Sample_", meta$ID)

bulk <- bulk[, meta$ID]
# bulk <- as.matrix(bulk)

fwrite(bulk, "/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz", row.names = T)

# bulk <- as.data.frame(fread("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz"))
# colnames(bulk)[1] <- "gene"
# row.names(bulk) <- bulk$gene
# bulk$gene <- NULL
# bulk <- as.matrix(bulk)
# bulk[1:10, 1:10]

# Load the Proportion results
prop_res <- as.data.frame(fread("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/res/Pancreas/proportion/scaden_all.tsv"))
prop_res$V1 <- paste0("Sample_", prop_res$V1)
row.names(prop_res) <- prop_res$V1
prop_res$V1 <- NULL

# Load the scData
scRNA <- zellkonverter::readH5AD("/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/diabetes_training.h5ad", reader = "R") # 124583  21487
assay(scRNA, "counts") <- assay(scRNA, "X")
scRNA <- scRNA[, scRNA$disease_state == "Control"]
scRNA_deconv <- scRNA[, scRNA$donor_id == "HPAP101"]
zellkonverter::writeH5AD(scRNA_deconv, "/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/pancreas_training.h5ad", X_name = "counts")
table(scRNA_deconv$CellType)

# library(CausalCellInfer)
# res <- infer_CTS_profile("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz",
#                          "/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/pancreas_training.h5ad",
#                          "/mnt/md0/yujia/project/github_package/demo/res/Pancreas/proportion/scaden_all.tsv",
#                          alpha=0.5, beta=0.1, tao_k=0.01, max.iter=1000,
#                          do_cpm = T, model_tracker = F, verbose = T,
#                          model_name = "model_1", Normalize = F, preprocess = "none", pos = F)

# bulk_dat <- "/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz"
# sc_dat <- "/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/pancreas_training.h5ad"
# prop <- "/mnt/md0/yujia/project/github_package/demo/res/Pancreas/proportion/scaden_all.tsv"

egm <- create_ENIGMA(bulk = bulk, ref = scRNA_deconv, ref_type = "single_cell")
egm <- batch_correct(egm, varname_cell_type = "CellType", ncores = 40)
prop_res <- prop_res[, colnames(egm@ref)]
egm@result_cell_proportion <- as.matrix(prop_res[colnames(bulk), ])

egm <- ENIGMA_L2_max_norm(egm, 
                          alpha=0.5, beta=0.1, tao_k=0.01, max.iter=1000,
                          do_cpm = T, model_tracker = F, verbose = T,
                          model_name = "model_1", Normalize = F, preprocess = "none", pos = F)
tmp_matrix <- assay(egm@result_CSE, "counts")
colnames(tmp_matrix) <- colData(egm@result_CSE)$label
tmp_meta <- colData(egm@result_CSE) 
patient_ids <- unique(tmp_meta$sample)
arr2 <- abind(lapply(split(seq_len(ncol(tmp_matrix)), 
            (seq_len(ncol(tmp_matrix))-1) %/% length(colnames(egm@ref)) + 1), function(x) tmp_matrix[, x]), along = 3)
tmp_res <- arr2

cell_labels <- colnames(egm@ref)
for (cell_index in 1:length(cell_labels)){
    # cell_index <- 1
    tmp_cell_exp <- tmp_res[,cell_index,]
    colnames(tmp_cell_exp) <- patient_ids
    tmp_cell_exp <- as.data.frame(tmp_cell_exp)

    fwrite(tmp_cell_exp, paste0("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/res/Pancreas/CTS/", stringr::str_replace(cell_labels[cell_index], " ", "_"), ".csv.gz"), 
            compress = "gzip", row.names = TRUE, col.names = TRUE)
}

library(data.table)
library(xgboost)
library(stringr)
library(zellkonverter)
library(scater)

source("/mnt/md0/yujia/project/github_package/CausalCellInfer/bin/cell_marker_identification/get_cellmarkers.R")

# Read the sc data
# scRNA <- fread("/home/yinly/Cell_deconvolution/01_Cellmarker_identification/data/T1D_T2D.csv") # dim(scRNA): 222077  21487
scRNA <- zellkonverter::readH5AD("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/sc_dat/Pancreas/T1D_T2D_public.h5ad", reader = "R") # 222077  21487
assay(scRNA, "counts") <- assay(scRNA,"X")

# Log-Normalization
scRNA <- logNormCounts(scRNA)

unknown_index_overall <- which(scRNA$cell_type=='unknown')
length(unknown_index_overall)

# Remove cell type with unknown annotation
scRNA <- scRNA[, scRNA$cell_type!='unknown']

# Select the donorID that will not be used for future deconvolution training.
donor_ids <- unique(scRNA$donor_id)
# external ids for T1D,T2D, control
#                  T1D       T1D       T2D       T2D       Control   Control   AAB       AAB
external_ids <- c("HPAP084","HPAP087","HPAP106","HPAP109","HPAP104","HPAP105","HPAP092","HPAP107")
selected_ids <- setdiff(donor_ids,external_ids)
scRNA <- scRNA[, scRNA$donor_id %in% selected_ids]

# Revised by yinly: August 21, 2023
# Description: we only include T2D and control for the analysis
T2D_index <- which(scRNA$disease_state=="T2D")
ctrl_index <- which(scRNA$disease_state=="Control")
selected_index <- union(T2D_index,ctrl_index)
scRNA_selected <- scRNA[, selected_index]

# classify the data into 2 different envirs based on diagnosis status
female_index <- which(scRNA_selected$sex=="female")
male_index <- which(scRNA_selected$sex=="male")

celltype_all <- as.vector(unlist(scRNA_selected$cell_type))
unique(celltype_all)

female <- scRNA_selected[, female_index]
male <- scRNA_selected[, male_index]

# use xgboost to build a prediction model and extract variable importance for T2D
data_fir <- t(as.matrix(assay(female, "logcounts")))
class(data_fir) <- "numeric"
label_fir <- as.vector(unlist(female$cell_type))

# use xgboost to build a prediction model and extract variable importance for controls
data_sec <- t(as.matrix(assay(male, "logcounts")))
class(data_sec) <- "numeric"
label_sec <- as.vector(unlist(male$cell_type))
rm(scRNA)

data_1 <- get_cellmarkers(label_fir, data_fir)
data_2 <- get_cellmarkers(label_sec, data_sec)

cellmarkers_common <- list()
cellmarkers_all <- list()

for (i in names(data_1[[2]])){

    cellmakers_fir <- as.data.frame(data_1[[2]][i][[1]])[, 1]
    cellmakers_sec <- as.data.frame(data_2[[2]][i][[1]])[, 1]
    celltype_common <- intersect(cellmakers_fir,cellmakers_sec)
    celltype_all <- union(cellmakers_fir,cellmakers_sec)
    cellmarkers_common[[i]] <- celltype_common
    cellmarkers_all[[i]] <- celltype_all
    # break
}

cellmarkers_all <- unique(unlist(cellmarkers_all))
cellmarkers_all.df <- data.frame("gene" = cellmarkers_all)
dim(cellmarkers_all.df)

cellmarkers_common <- unique(unlist(cellmarkers_common))
cellmarkers_common.df <- data.frame("gene" = cellmarkers_common)
dim(cellmarkers_common.df)

fwrite(cellmarkers_all.df, "/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/marker_genes/Pancreas/all_genes.csv")
fwrite(cellmarkers_common.df, "/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/marker_genes/Pancreas/common_genes.csv")

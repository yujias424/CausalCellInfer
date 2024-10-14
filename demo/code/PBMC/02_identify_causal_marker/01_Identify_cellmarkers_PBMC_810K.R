library(data.table)
library(xgboost)
library(stringr)

source("/mnt/md0/yujia/project/github_package/CausalCellInfer/bin/cell_marker_identification/get_cellmarkers.R")

scRNA <- fread("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/sc_dat/pbmc810k_normalizedX.csv.gz")
cell_info <- fread("/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/sc_dat/pbmc810k_meta.csv")
celltype <- cell_info$cell_type
IID <-cell_info$IID

scRNA_all <- cbind(scRNA,celltype,IID)

# June 4, 2024
# Description: remove Dendritic Cells to keep it consistent with Newman dataset!
# known_index <- which(celltype!="Megakaryocytes")
mege_index <- which(celltype=="Megakaryocytes")
dc_index <- which(celltype=="Dendritic Cells")
unknown_index <- union(mege_index,dc_index)
known_index <- setdiff(seq(1,length(celltype)),unknown_index)
scRNA_selected <- scRNA_all[known_index,]

feat_len <- ncol(scRNA_all)- 2

donor10K_index <- which(scRNA_selected$IID=="pbmc10k")
donor10K <- scRNA_selected[donor10K_index,]

donor8K_index <- which(scRNA_selected$IID=="pbmc8k")
donor8K <- scRNA_selected[donor8K_index,]

data_fir <- subset(donor10K,select=(1:feat_len))
data_fir <- as.matrix(data_fir)
class(data_fir) <- "numeric"
label_fir <- as.vector(unlist(donor10K$celltype))

data_sec <- subset(donor8K,select=(1:feat_len))
data_sec <- as.matrix(data_sec)
class(data_sec) <- "numeric"
label_sec <- as.vector(unlist(donor8K$celltype))

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

fwrite(cellmarkers_all.df, "/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/marker_genes/all_genes.csv")
fwrite(cellmarkers_common.df, "/mnt/md0/yujia/project/github_package/CausalCellInfer/demo/dat/marker_genes/common_genes.csv")

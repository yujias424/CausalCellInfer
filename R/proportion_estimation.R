scaden_py <- NULL
pd <- NULL

#' Load the scaden_py module for deconvolution
#'
#' The .onLoad code is used to load the scaden_py module for deconvolution,
#'
#' @param sc_dat the reference single cell data.
#' @param bulk_dat the bulk data to be deconvoluted.
#'
#' @return a result matrix with 3 columns, i.e., Genes, zMin, Assoc
#'
#' @importFrom reticulate import_from_path import
#' 
#' @export
.onLoad <- function(libname, pkgname) {
    # use superassignment to upandasate global reference to scipy
    path <- system.file("python", package = "CausalCellInfer")
    scaden_py <<- import_from_path("scaden_py", path = path, delay_load = TRUE)
    pd <<- import("pandas", as = "pd", delay_load = TRUE)
}

#' Obtain the proportion of the cell types of a target 
#'
#' The following code are used to identify cell type specific causal genes for the target trait
#'
#' @param sc_dat the reference single cell data.
#' @param bulk_dat the bulk data to be deconvoluted (Gene at row and Sample at column)
#' @param marker_genes marker genes for deconvolution (A csv file with a column named "gene")
#' @param ... arguments from the scaden deconvolution.
#'
#' @return a result matrix with 3 columns, i.e., Genes, zMin, Assoc
#'
#' @import reticulate
#' 
#' @export
get_proportion <- function(sc_dat, bulk_dat, marker_genes, sep, ...){

    # Read the bulk data
    bulk_dat <- pd$read_csv(bulk_dat, index_col = "Unnamed: 0", sep=sep)$T # We expect gene at row and sample at column
    marker_genes <- read.csv(marker_genes)$gene # We expect a csv file with a specific gene column

    pred_Pancreas_common <- scaden_py$ScadenDeconvolution(sc_dat, bulk_dat, markergene = marker_genes, ...)

    pred_Pancreas_common
} 


# get_proportion <- function(sc_dat, bulk_dat, marker_genes){

#     # Read the bulk data
#     bulk_dat <- pandas$read_csv(bulk_dat, index_col = "Unnamed: 0", sep="\t")$T

#     pred_Pancreas_common <- scaden_py$ScadenDeconvolution(sc_dat, bulk_dat,
#                                                           generate_sim_method = "tape", add_noise = FALSE, cut_variance= TRUE, cell_type = "CellType",
#                                                           markergene = marker_genes)

#     pred_Pancreas_common
# }


# bulk_dat <- "/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz"
# use_condaenv("~/miniconda3/envs/causalcellinfer")
# pandas <- import("pandas", as = "pandas")
# bulk_dat_aa <- pandas$read_csv(bulk_dat,   sep=",") # index_col = "Unnamed: 0",
# marker_genes <- read.csv("/mnt/md0/yujia/project/github_package/demo/dat/marker_genes/Pancreas/common_genes.csv")$gene

# get_proportion("/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/pancreas_training.h5ad", "/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz", "/mnt/md0/yujia/project/github_package/demo/dat/marker_genes/Pancreas/common_genes.csv",
#                 generate_sim_method = "tape", add_noise = FALSE, cut_variance= TRUE, cell_type = "CellType")

# library(reticulate)
# library(CausalCellInfer)
# use_condaenv("~/miniconda3/envs/causalcellinfer")
# pandas <- import("pandas", as = "pandas")
# scaden_py <- import_from_path("scaden_py", path = "/mnt/md0/yujia/project/github_package/CausalCellInfer/inst/python/")

# # bulk_dat <- pandas$read_csv("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_example.txt.gz", sep = "\t", index_col = "Unnamed: 0")$T
# bulk_dat <- pandas$read_csv("/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_bulk.txt.gz", sep = ",", index_col = "Unnamed: 0")$T
# common_genes <- read.csv("/mnt/md0/yujia/project/github_package/demo/dat/marker_genes/Pancreas/common_genes.csv")$gene

# pred_Pancreas_common <- scaden_py$ScadenDeconvolution("/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/pancreas_training.h5ad", bulk_dat,
#                                            generate_sim_method = "tape", add_noise = FALSE, cut_variance= TRUE, cell_type = "CellType",
#                                            markergene = common_genes)

# library(CausalCellInfer)
# reticulate::use_condaenv("~/miniconda3/envs/causalcellinfer")

# marker_genes <- read.csv("/mnt/md0/yujia/project/github_package/demo/dat/marker_genes/Pancreas/common_genes.csv")$gene
# res <- get_proportion(sc_dat = "/mnt/md0/yujia/project/github_package/demo/dat/sc_dat/Pancreas/diabetes_training.h5ad",
#                       bulk_dat = "/mnt/md0/yujia/project/github_package/demo/dat/bulk/Pancreas/GSE50244_example.txt.gz",
#                       marker_genes = marker_genes,
#                       generate_sim_method = "tape", add_noise = FALSE, cut_variance= TRUE, cell_type = "CellType")

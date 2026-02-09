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
    bulk_dat <- pd$read_csv(bulk_dat, index_col = as.integer(0), sep=sep)$T # We expect gene at row and sample at column, noticed that we need as.integer here, otherwise the reticulate will treat it as 0.0
    marker_genes <- read.csv(marker_genes)$gene # We expect a csv file with a specific gene column

    pred_Pancreas_common <- scaden_py$ScadenDeconvolution(sc_dat, bulk_dat, markergene = marker_genes, ...)

    pred_Pancreas_common
} 

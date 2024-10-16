#' Infer the CTS profile using ENIGMA
#'
#' The following code is to use the ENIGMA to infer the CTS profile for different cell types.
#'
#' @param bulk_dat bulk RNA-seq dataset, we expect the bulk data index should match the propotion data
#' @param sc_dat single cell RNA-seq dataset.
#' @param prop estimated proportion of the bulk RNA-seq dataset.
#' 
#' @return a list containing the estimated CTS profile
#'
#' @import ENIGMA
#' @import scater
#' @import abind
#' @importFrom zellkonverter readH5AD
#' @importFrom data.table fread
#' 
#' @export
infer_CTS_profile <- function(bulk_dat, sc_dat, prop, ...){

    # Read the bulk data (We expect sample at column, gene at row)
    bulk <- as.data.frame(fread(bulk_dat))
    colnames(bulk)[1] <- "gene"
    row.names(bulk) <- bulk$gene
    bulk$gene <- NULL
    bulk <- as.matrix(bulk)

    # Read the proportion data
    prop_res <- as.data.frame(fread(prop))
    prop_res$V1 <- paste0("Sample_", prop_res$V1)
    row.names(prop_res) <- prop_res$V1
    prop_res$V1 <- NULL

    # Load the scData
    scRNA <- readH5AD(sc_dat, reader = "R") # We assume the counts assay is exists.
    
    # Create the ENIGMA object
    egm <- create_ENIGMA(bulk = bulk, ref = scRNA, ref_type = "single_cell")
    egm <- batch_correct(egm, varname_cell_type = "CellType", ncores = 40)
    prop_res <- prop_res[, colnames(egm@ref)]
    # egm@result_cell_proportion <- as.matrix(prop_res[colnames(bulk), ])
    egm@result_cell_proportion <- as.matrix(prop_res)

    # ENIGMA inference
    egm <- ENIGMA_L2_max_norm(egm, ...)

    # Summarize the estimated CTS profile.
    tmp_matrix <- assay(egm@result_CSE, "counts")
    colnames(tmp_matrix) <- colData(egm@result_CSE)$label
    tmp_meta <- colData(egm@result_CSE) 
    patient_ids <- unique(tmp_meta$sample)
    arr2 <- abind(lapply(split(seq_len(ncol(tmp_matrix)), 
                (seq_len(ncol(tmp_matrix))-1) %/% length(colnames(egm@ref)) + 1), function(x) tmp_matrix[, x]), along = 3)
    tmp_res <- arr2

    cell_labels <- colnames(egm@ref)
    res_list <- list()

    for (cell_index in 1:length(cell_labels)){
        tmp_cell_exp <- tmp_res[,cell_index,]
        colnames(tmp_cell_exp) <- patient_ids
        tmp_cell_exp <- as.data.frame(tmp_cell_exp)

        res_list[[cell_labels[cell_index]]] <- tmp_cell_exp
    }

    res_list
}
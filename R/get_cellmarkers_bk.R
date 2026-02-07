# #' Get the cell markers for cell deconvolution using ICP approach
# #'
# #' The following code are used to identify environment specific cellmarkers.
# #'
# #' @param label the corresponding cell type labels for the input data
# #' @param data input normalized scRNA-seq data(count data would also be OK), numeric data type
# #'
# #' @return A list object containing all target cell types and their corresponding cell markers.
# #'
# #' @importFrom xgboost xgb.DMatrix xgb.train xgb.importance
# #' @import data.table
# #' 
# #' @export
# get_cellmarkers <- function(label,data){
#     celltypes <- unique(label)
#     marker_list <- list()
#     for (i in 1:length(celltypes)){
#         celltype <- celltypes[i]
#         label_binary <- rep(0,length(label))
#         celltype_index <- which(label==celltype)
#         cat("There are",length(celltype_index),celltype,"cells in the current dataset. \n")
#         label_binary[celltype_index] <- 1
#         dtrain<- xgb.DMatrix(data = data, label = label_binary)
#         dtrain_model <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 20, nrounds = 10, objective = "binary:logistic")
#         dtrain_rank <- xgb.importance(model = dtrain_model)
#         cat("There are",nrow(dtrain_rank),"identified cellmarkers for the current dataset.\n")
#         marker_list[[celltype]] <- dtrain_rank  
#         rm(dtrain, dtrain_model, dtrain_rank)
#         gc()
#     }
#     result_list <- list(celltypes=celltypes,marker_list=marker_list)
#     return(result_list)
# }


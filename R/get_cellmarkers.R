#' Re-rank genes based on zMin
#' 
#' merge_and_rank_genes is designed to re-rank genes based on zMin
#' 
#' @param df1 dataframe with three columns: genes, assoc and zMin to store the identifed caual genes from environment 1 
#' @param df2 dataframe with three columns: genes, assoc and zMin to store the identifed caual genes from environment 2
#' 
#' @return final_df re-rank the merged df in descending order(zMin), with zMin indicates variable importance
#' 
#' @export 
merge_and_rank_genes <- function(df1, df2) {

    # library(data.table)
    # library(xgboost)
    # library(GeneralisedCovarianceMeasure)
    # library(pcalg)
    # source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/I-GCM.R") # This will be packed into a R pacakge
    # source("/mnt/home/yinly/projects/Cell_deconvolution/bin/src/PCSelect_Parallel_Update.R")
    #**************************************************************************************************************************
    #' @ Description: merge_and_rank_genes is designed to re-rank genes based on zMin
    #' @para df1: dataframe with three columns: genes, assoc and zMin to store the identifed caual genes from environment 1 
    #' @para df2: dataframe with three columns: genes, assoc and zMin to store the identifed caual genes from environment 2
    #' @return final_df: re-rank the merged df in descending order(zMin), with zMin indicates variable importance
    #**************************************************************************************************************************
    # Ranked genes by zMin

    # Ensure required columns exist
    required_cols <- c("genes", "assoc", "zMin")
    stopifnot(all(required_cols %in% colnames(df1)),
                all(required_cols %in% colnames(df2)))
    
    # Merge by gene
    merged <- merge(df1, df2, by = "genes", all = TRUE, suffixes = c(".1", ".2"))
    
    # Compute average zMin (NA-safe)
    merged$zMin <- rowMeans(cbind(as.numeric(merged$zMin.1), as.numeric(merged$zMin.2)), na.rm = TRUE)
    
    # Combine assoc as logical OR (TRUE if TRUE in either)
    merged$assoc <- (merged$assoc.1 %in% TRUE) | (merged$assoc.2 %in% TRUE)
    
    # Create final output
    final_df <- merged[, c("genes", "assoc", "zMin")]
    
    # Sort by zMin decreasing
    final_df <- final_df[order(-final_df$zMin), ]
    
    rownames(final_df) <- NULL
    return(final_df)
}


#' Get cellmarkers using a combination of PC-simple and I-GCM
#' 
#' get_cellmarkers is designed to automatically identify cellmarkers using a combination of PC-simple and I-GCM
#' 
#' @param scRNA the scRNA-seq data including both environment
#' @param celltypes_all the corresponding cell type lable for scRNA
#' @param envir_binary the binary environment variable
#' @param data_fir scRNA-seq data from the 1st environment
#' @param label_fir celltype labels corresponding to data_fir
#' @param data_sec scRNA-seq data from the 2nd environment
#' @param label_sec celltype labels corresponding to data_sec
#' @param celltypes a vector contain all unique cell types
#' @param sample_num overall sample size for scRNA-seq from both environments
#' @param alpha the pvalue cutoff used to detect causal variables in PC-simple algorithm, the default value is 0.05
#' @param change_alpha the pvalue cutoff used to detect GCM distance change significance, the default value is 0.05
#' @param verbose indicates whether outputS the analysis details, either TRUE or FALSE
#' 
#' @return res_list contains results list for each cell types. There are 6 variables in each sublist, i.e., celltype, res_causal,candidate_causal_genes,IGCM_mat,causal_set_index,causal_genes_ICP, which separately specify the celltype, candidate variable matrix with their corresponding zMin, candidate causal genes, the IGCM detection results matrix, causal variable set index(the last index of the variable),the identified causal genes using ICP 
#' 
#' @import data.table,
#' @import xgboost
#' @import GeneralisedCovarianceMeasure
#' @import pcalg
#' 
#' @export 

get_cellmarkers <- function(scRNA,celltypes_all,envir_binary,data_fir,label_fir,data_sec,label_sec,celltypes,sample_num,alpha=0.05,change_alpha=0.05,verbose=FALSE){

    #*************************************************************************************************
    #' @ Description: get_cellmarkers is designed to automatically identify cellmarkers using a 
    # combination of PC-simple and I-GCM
    #' @para scRNA: the scRNA-seq data including both environment
    #' @para celltypes_all: the corresponding cell type lable for scRNA
    #' @para envir_binary: the binary environment variable
    #' @para data_fir: scRNA-seq data from the 1st environment
    #' @para label_fir: celltype labels corresponding to data_fir
    #' @para data_sec: scRNA-seq data from the 2nd environment
    #' @para label_sec: celltype labels corresponding to data_sec
    #' @para celltypes: a vector contain all unique cell types
    #' @para sample_num: overall sample size for scRNA-seq from both environments
    #' @para alpha: the pvalue cutoff used to detect causal variables in PC-simple algorithm, the default value is 0.05
    #' @para change_alpha: the pvalue cutoff used to detect GCM distance change significance, the default value is 0.05
    #' @para verbose: indicates whether outputS the analysis details, either TRUE or FALSE
    #' @return res_list: contains results list for each cell types. There are 6 variables in each sublist, 
    #' i.e., celltype, res_causal,candidate_causal_genes,IGCM_mat,causal_set_index,causal_genes_ICP, which separately specify
    #' the celltype, candidate variable matrix with their corresponding zMin, candidate causal genes, the IGCM detection results matrix, causal variable set index(the last index of the variable),the identified causal genes using ICP 
    #*************************************************************************************************

    res_list <- list()
    for(i in 1:length(celltypes)){
        celltype <- celltypes[i]
        if(verbose){
            cat("This is the identification of cell marker for",celltype,"\n")
        }
        # train xgboost for the 1st dataset
        label_fir_binary <- rep(0,length(label_fir))
        celltype_index <- which(label_fir==celltype)
        if(verbose){
            cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
        }   
        label_fir_binary[celltype_index] <- 1
        dtrain_fir <- xgb.DMatrix(data = data_fir, label = label_fir_binary)
        dtrain_fir_model <- xgb.train(data = dtrain_fir, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
        dtrain_fir_rank <- xgb.importance(model = dtrain_fir_model)
        if(verbose){
            cat("There are",nrow(dtrain_fir_rank),"identified cellmarkers in the 1st dataset.\n")
        }
        assoc_feats <- as.vector(unlist(dtrain_fir_rank$Feature))
        X_assoc <- subset(data_fir,select=assoc_feats)
        # employ PCSimple to identify the causal cellmarkers
        dfnew = data.frame(outcome=label_fir_binary, X_assoc,check.names = FALSE) # binarized celltype labels
        library(coop)
        expr_corMat = pcor(dfnew[,-1])
        #now add back the correlation between the outcome and each feature
        #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
        cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)}  )
        precompute_corMat = cbind(cor_with_outcome, expr_corMat)
        precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat)
        n <- nrow (dfnew)
        V <- colnames(dfnew) # labels aka node names
        ## estimate local causal network structure around the response variable
        #library(parallel)
        pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
                                            dm=dfnew[,-1],
                                            method = c("parallel"),
                                            mem.efficient = FALSE,
                                            num_workers = 10,
                                            # alpha=0.001,
                                            alpha=alpha,
                                            corMat = precompute_corMat,
                                            max_ord=3,
                                            corMethod = "standard",
                                            verbose = TRUE, directed = TRUE)
        # The following codes are used to summarize causally relevant genes
        assoc <- pcSimple.fit$G
        genes <- names(pcSimple.fit$G)
        # genes <- gsub("HLA.","HLA-",genes) # There are special cases where the "-" in gene symbol was replaced by ".
        zMin <- pcSimple.fit$zMin
        res <- cbind(genes,assoc,zMin)
        res_causal <- res[which(assoc==TRUE),]
        if(!is.null(dim(res_causal))){
            causal_genes <- res_causal[,1]
        }else{
            causal_genes <- res_causal[1]
        }
        
        # Train xgboost for the 2nd dataset
        label_sec_binary <- rep(0,length(label_sec))
        celltype_index <- which(label_sec==celltype)
        if(verbose){
            cat("There are",length(celltype_index),celltype,"cells in the 1st dataset. \n")
        }       
        label_sec_binary[celltype_index] <- 1
        dtrain_sec <- xgb.DMatrix(data = data_sec, label = label_sec_binary)
        dtrain_sec_model <- xgb.train(data = dtrain_sec, max.depth = 2, eta = 1, nthread = 10, nrounds = 10, objective = "binary:logistic")
        dtrain_sec_rank <- xgb.importance(model = dtrain_sec_model)
        if(verbose){
            cat("There are",nrow(dtrain_sec_rank),"identified cellmarkers in the 2nd dataset.\n")
        }
        assoc_feats_sec <- as.vector(unlist(dtrain_sec_rank$Feature))
        X_assoc_sec <- subset(data_sec,select=assoc_feats_sec)
        # employ PCSimple to identify the causal cellmarkers
        dfnew = data.frame(outcome=label_sec_binary, X_assoc_sec,check.names=FALSE) # binarized celltype labels
        # library(coop)
        expr_corMat = pcor(dfnew[,-1])
        #now add back the correlation between the outcome and each feature
        #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
        cor_with_outcome = apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)}  )
        precompute_corMat = cbind(cor_with_outcome, expr_corMat)
        precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat)
        n <- nrow (dfnew)
        V <- colnames(dfnew) # labels aka node names
        ## estimate local causal network structure around the response variable
        #library(parallel)
        pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
                                            dm=dfnew[,-1],
                                            method = c("parallel"),
                                            mem.efficient = FALSE,
                                            num_workers = 10,
                                            # alpha=0.001,
                                            alpha=alpha,
                                            corMat = precompute_corMat,
                                            max_ord=3,
                                            corMethod = "standard",
                                            verbose = TRUE, directed = TRUE)
        # The following codes are used to summarize causally relevant genes
        assoc <- pcSimple.fit$G
        genes <- names(pcSimple.fit$G)
        # genes <- gsub("HLA.","HLA-",genes) # There are special cases where the "-" in gene symbol was replaced by "."
        # genes <- gsub("MT.","MT-",genes) # There are special cases where the "-" in gene symbol was replaced by "."
        zMin <- pcSimple.fit$zMin
        res_sec <- cbind(genes,assoc,zMin)
        res_causal_sec <- res_sec[which(assoc==TRUE),]
        if(!is.null(dim(res_causal_sec))){
            causal_genes_sec <- res_causal_sec[,1]
        }else{
            causal_genes_sec <- res_causal_sec[1]
        }
        
        res_causal_merged <- merge_and_rank_genes(res_causal,res_causal_sec)
        if(!is.null(dim(res_causal_merged))){
            causal_genes_final <- res_causal_merged[,1]
        }else{
            causal_genes_final <- res_causal_merged[1]
        }
                
        # Employing PC
        # Y_binary <- label_binary
        Y_binary <- rep(0,length(celltypes_all))
        celltype_index <- which(celltypes_all==celltype)
        Y_binary[celltype_index] <- 1 # celltype label: one vs remaining all

        X_assoc <- subset(scRNA,select=causal_genes_final)
        initial_assoc_len <- ncol(X_assoc)
        # alpha <-0.05
        res_mat <- matrix(1,nrow=ncol(X_assoc),ncol=3)
        colnames(res_mat)<- c("variables","stat","p")
    
        resid_mat <- NULL
        for(j in initial_assoc_len:1){
            X<-envir_binary
            present_index_set <- 1:j
            Z <- subset(X_assoc,select=present_index_set)
            gcm_obj <- gcm_test(X=X,Y=Y_binary,Z=Z,regr.par=list(max_depth=4,max_delta_step=6), regr.method = "xgboost")
            #gcm_obj <- gcm_test(X=X,Y=Y_binary,Z=Z,regr.par=list(max_depth=3,max_delta_step=10), regr.method = "xgboost")
            gcm_obj_stat <- gcm_obj$test.statistic
            gcm_obj_p <- gcm_obj$p.value
            gcm_R <- gcm_obj$R
            resid_mat <- cbind(resid_mat,gcm_R)
            #cat("The statistic and pvalue for independence test are respectively",gcm_obj_stat,"and",gcm_obj_p,"\n")
            res_mat[(initial_assoc_len-j+1),1] <- j
            res_mat[(initial_assoc_len-j+1),2] <- gcm_obj_stat
            res_mat[(initial_assoc_len-j+1),3] <- gcm_obj_p
        }
        IGCM_mat <- matrix(0,nrow=nrow(res_mat)-1,ncol=5)
        colnames(IGCM_mat) <- c("base_set","reduced_set","stat_change","stat_change_Z","stat_change_P")
        causal_set_index <- NULL
        # change_alpha <- 0.05
        for(k in 1:(nrow(res_mat)-1)){
            base_gcm_stat <- res_mat[k,2]
            next_gcm_stat <- res_mat[(k+1),2]
            R1 <- resid_mat[,k]
            R2 <- resid_mat[,(k+1)]
            gcm_stat_stat_cov <- stat_change_detection(R1,R2,sample_num) # for Tns that follows standard normal distribution, the correlation also equals the covariance
            gcm_stat_change <- abs(next_gcm_stat) - abs(base_gcm_stat)
            gcm_stat_change_var <- 2-4/pi*(gcm_stat_stat_cov*asin(gcm_stat_stat_cov)+sqrt(1-gcm_stat_stat_cov^2))
            if(abs(gcm_stat_change)>0){
                gcm_stat_change_Z <- gcm_stat_change/sqrt(gcm_stat_change_var) # we directly use the calculated change!
                gcm_stat_change_P <- pnorm(gcm_stat_change_Z,lower.tail = FALSE) # Here, we use one tail test
                IGCM_mat[k,1] <- nrow(res_mat)-k+1
                IGCM_mat[k,2] <- nrow(res_mat)-k
                IGCM_mat[k,3] <- gcm_stat_change
                IGCM_mat[k,4] <- gcm_stat_change_Z
                IGCM_mat[k,5] <- gcm_stat_change_P
                if(verbose){
                    cat("The stats change between the first",nrow(res_mat)-k+1,"and",nrow(res_mat)-k,"variables is",gcm_stat_change,"with a significance of",gcm_stat_change_P,"\n")                
                }
                if(gcm_stat_change_P<=change_alpha){
                    if(verbose){
                        cat("The first",nrow(res_mat)-k+1,"variables is the potential causal variable set.\n")
                    }                 
                    causal_index <- nrow(res_mat)-k+1
                    causal_set_index <- c(causal_set_index,causal_index)
                }
            }else{
                    if(verbose){
                        cat("The stats change between the first",nrow(res_mat)-k+1,"and",nrow(res_mat)-k,"variables is",gcm_stat_change,"\n")
                    }                 
            }
        }       
        marker_gene_index <- causal_set_index[1]
        markers <- causal_genes_final[1:marker_gene_index]
        # outputfile <- paste0("/mnt/home/yinly/projects/Cell_deconvolution/res/06_Benchmark/01_Cellmarker_identification/PBMC/6810KABC_Batch_PCSimple_ICP_New/6810KABC_cellmarkers_",celltype,".RData")
        # save(celltype,dtrain_fir_rank,dtrain_sec_rank,res_causal,res_causal_sec,res_causal_merged,IGCM_mat,causal_genes_final,causal_set_index,file=outputfile)
        cat("The identification of cellmarkers for",celltype ,"is completed!\n\n")
        causal_genes_ICP_index <- causal_set_index[1]
        if(length(causal_set_index)>=1){
            causal_genes_ICP <- causal_genes_final[1:causal_genes_ICP_index]
        }else{
            causal_genes_ICP <- NULL
        }
        res_celltype_list <- list(
            res_causal = res_causal_merged,
            candidate_causal_genes = causal_genes_final,
            IGCM_mat = IGCM_mat,
            causal_set_index = causal_set_index,
            causal_genes_ICP = causal_genes_ICP
        )
        res_list[[celltype]] <- res_celltype_list  
    }
    return(res_list)  
}

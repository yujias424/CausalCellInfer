#' Calculate the covariance between the calculated GCMs from two different GCM tests
#' 
#' The following code is used to calculate the covariance between the calculated GCMs from two different GCM tests.
#' 
#' @param R1 derived residual product for the 1st test
#' @param R2 derived residual product for the 2nd test
#' @param N sample size
#' 
#' @return calculated covariance
#' 
#' @import data.table
#' @import coop
#' @import parallel
#' @import pcalg
#' @import GeneralisedCovarianceMeasure
#' 
#' @export 
gcm_cov_detection <- function(R1, R2, N){
    
    # **********************************************************************************************************************************************
    # @ Description: calculate the covariance between the calculated GCMs from two different GCM tests
    # @ R1: derived residual product for the 1st test
    # @ R2: derived residual product for the 2nd test
    # @ N: sample size
    # @ tn_cov: calculated covariance
    # library(data.table)
    # library(coop)
    # library(parallel)
    # library(pcalg)
    # library(GeneralisedCovarianceMeasure)
    # source("./src/PCSelect_Parallel_Update.R")
    # **********************************************************************************************************************************************

    se_R1 <- sqrt(var(R1))/sqrt(N)
    se_R2 <- sqrt(var(R2))/sqrt(N)
    tn_cov <- 1/(se_R1*se_R2) * cov(R1,R2)/N
    return(tn_cov)
}

#' I-GCM is designed to accurately identify direct causal variable set for the target variable
#' 
#' The invariantGCM function is to accurately identify direct causal variable set for the target variable
#' 
#' @param X_assoc all input variables excluding the environment and target variable, with each row represents one subject
#' @param Y target variable vector for the same set of subjects
#' @param envir environment variable vector for the same set of subjects
#' @param num_workers number of threads used for feature selection
#' @param change_alpha p-value cutoff used to detect significant GCM change, default alpha=0.05
#' @param change_gcm absolute gcm change threshold for significant change, default change_gcm=0
#' @param verbose if verbose=TRUE, the details of the analysis will be print out. The default value is false
#'
#' @return res a list stored all related results, which includes candidate_features, resid_mat and all_candidate_causal_set; candidate_features: the preselected candidate causal variables from PC-simple ranked by zscore in descending order; resid_mat: the residual matrix with each column stores the residuals for all subjects in each performed test; all_candidate_causal_set: number of causal variables(The top ranked variables from candidate_features) for each identified 
#'
#' @import pcalg
#' @importFrom coop pcor
#' 
#' @export
InvariantGCM <- function(X_assoc, Y, envir, change_alpha=0.05,
                         change_gcm=0, num_workers, verbose=FALSE){

    #**********************************************************************************************************************************************
    # @ Description: I-GCM is designed to accurately identify direct causal variable set for the target variable
    # Notes: Make sure each row
    # @ X_assoc: all input variables excluding the environment and target variable, with each row represents one subject
    # @ Y: target variable vector for the same set of subjects
    # @ envir: environment variable vector for the same set of subjects
    # @ num_workers: number of threads used for feature selection
    # @ change_alpha: p-value cutoff used to detect significant GCM change, default alpha=0.05
    # @ change_gcm: absolute gcm change threshold for significant change, default change_gcm=0
    # @ verbose: if verbose=TRUE, the details of the analysis will be print out. The default value is false
    # @ res: a list stored all related results, which includes candidate_features,resid_mat and all_candidate_causal_set 
    # candidate_features: the preselected candidate causal variables from PC-simple ranked by zscore in descending order
    # resid_mat: the residual matrix with each column stores the residuals for all subjects in each performed test
    # all_candidate_causal_set: number of causal variables(The top ranked variables from candidate_features) for each identified causal variable set  
    #**********************************************************************************************************************************************

    Y_outcome <- Y
    dfnew <- data.frame(outcome = Y_outcome, X_assoc)
    dfnew_rankNorm <- dfnew
    expr_corMat <- pcor(dfnew[,-1])
    cor_with_outcome <- apply(dfnew[,-1],  2 , function(col) {pcor(dfnew[,"outcome"], col)})
    precompute_corMat <- cbind(cor_with_outcome, expr_corMat)
    precompute_corMat <- rbind(c(1,cor_with_outcome), precompute_corMat)
    pcSimple.fit <- PCSelect_Parallel(y = dfnew[,1],
                                      dm = dfnew[,-1],
                                      method = c("parallel"),
                                      mem.efficient = FALSE,
                                      num_workers = num_workers,
                                      alpha = 0.05,
                                      corMat = precompute_corMat,
                                      max_ord = 3,
                                      corMethod = "standard",
                                      verbose = TRUE, directed = TRUE)
    feat_names <- names(pcSimple.fit$G)
    assoc <- as.vector(unlist(pcSimple.fit$G))
    zMin <- as.vector(unlist(pcSimple.fit$zMin))
    pcsimple_res <- cbind(feat_names,assoc,zMin)
    pcsimple_res_sort <- pcsimple_res[order(as.numeric(pcsimple_res[,3]),decreasing=TRUE),]
    assoc_feats <- pcsimple_res_sort[pcsimple_res_sort[,2]==TRUE,1]
    X_assoc <- X_assoc[,assoc_feats]
    if (verbose){
        cat("There are", length(assoc_feats), "variables remaining after PCSimple-based feature selection.\n")
    }
    initial_assoc_len <- ncol(X_assoc)
    resid_mat <- NULL # resid_mat is used to store the product of residuals for all subject
    for (i in initial_assoc_len:1){
        if (verbose){
            cat("This is the test with the first ", i, "variable.\n")
        } 
        present_index_set <- 1:i
        X <- envir
        Z <- X_assoc[,present_index_set]
        gcm_obj <- gcm_test(X = X, Y = Y_binary, Z = Z, regr.par = list(max_depth=4, max_delta_step=6), regr.method = "xgboost")
        gcm_obj_stat <- gcm_obj$test.statistic
        gcm_R <- gcm_obj$R
        gcm_obj_p <- gcm_obj$p.value
        if (verbose){
            cat("The statistic and pvalue for independence test are respectively", gcm_obj_stat, "and", gcm_obj_p, "\n")
        }
        resid_mat <- cbind(resid_mat,gcm_R)
        res_mat[(initial_assoc_len-i+1),1] <- i
        res_mat[(initial_assoc_len-i+1),2] <- gcm_obj_stat
        res_mat[(initial_assoc_len-i+1),3] <- gcm_obj_p
    }
    sample_num <- length(Y)
    causal_set_index <- NULL
    for (j in 1:(nrow(res_mat)-1)){
        base_gcm_stat <- res_mat[j, 2]
        next_gcm_stat <- res_mat[(j+1), 2]
        R1 <- resid_mat[, j]
        R2 <- resid_mat[, (j+1)]
        gcm_stat_stat_cov <- gcm_cov_detection(R1,R2,sample_num)
        gcm_stat_change <- abs(next_gcm_stat) - abs(base_gcm_stat) # Here we use the absolute deviation differences as the measure of significant change
        gcm_stat_change_var <- 2-4/pi*(gcm_stat_stat_cov*asin(gcm_stat_stat_cov)+sqrt(1-gcm_stat_stat_cov^2))
        if (abs(gcm_stat_change) > 0){
            gcm_stat_change_Z <- gcm_stat_change/sqrt(gcm_stat_change_var)
            gcm_stat_change_P <- pnorm(gcm_stat_change_Z,lower.tail = FALSE) # Here, we use one tail test
            if (verbose){
                cat("The stats change between the first", nrow(res_mat)-j+1, "and", nrow(res_mat)-j, "variables is",
                    gcm_stat_change, "with a significance of", gcm_stat_change_P, "\n")
            }
            if (gcm_stat_change_P <= change_alpha && abs(gcm_stat_change) > change_gcm){
                if (verbose){
                    cat("The first",nrow(res_mat)-j+1,"variables is the potential causal variable set.\n")
                }
                causal_index <- nrow(res_mat)-j+1
                causal_set_index <- c(causal_set_index,causal_index)
            }
        } else {
            cat("The stats change between the first", nrow(res_mat)-j+1, "and", 
                nrow(res_mat)-j, "variables is", gcm_stat_change, "please check your data.\n")
        }
    }
    res <- list(candidate_features = X_assoc, resid_mat = resid_mat, all_candidate_causal_set = causal_set_index)
    return(res)
}
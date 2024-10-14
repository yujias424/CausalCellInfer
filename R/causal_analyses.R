#' Get the causal genes using the proposed PC-simple algorithm
#'
#' The following code are used to identify cell type specific causal genes for the target trait
#'
#' @param data unconfounded cell type specific gene expression profiles for all genes
#' @param outcome the target trait
#' @param alpha the cutoff used to define causally relevant genes, we set alpha to 0.05 in this study, but you can adaptively adjust this cutoff
#' @param num_workers the number of cores used to calculate causally relevant genes, we set it to 10 in our study. Again you can adjust this to adapt to your device
#'
#' @return a result matrix with 3 columns, i.e., Genes, zMin, Assoc
#'
#' @importFrom WGCNA cor
#' 
#' @export
get_causal_genes <- function(data,outcome,alpha,num_workers){

    #******************************************************************************************************
    # The following code are used to identify cell type specific causal genes for the target trait
    # @funct get_causal_genes
    # @param data: unconfounded cell type specific gene expression profiles for all genes
    # @param outcome: the target trait
    # @param alpha: the cutoff used to define causally relevant genes, we set alpha to 0.05 in this study, but you can adaptively adjust this cutoff
    # @param num_workers: the number of cores used to calculate causally relevant genes, we set it to 10 in our study. Again you can adjust this to adapt to your device
    # @return result_mat: a result matrix with 3 columns, i.e., Genes, zMin, Assoc
    # Genes: lists the gene symbols of the causal genes
    # zMin: The minimal z-values when testing partial correlations between outcome and each gene. 
    # Assoc: Indicates whether the current gene is a causal gene based on predefined causal gene cutoff
    #******************************************************************************************************

    dfnew <- data.frame(outcome = outcome, data)
    expr_corMat <- cor(dfnew[,-1],use = "all.obs", method="pearson",nThreads = 20)
    cor_with_outcome <- apply(dfnew[,-1],  2 , function(col) {cor(dfnew[,"outcome"], col)})
    precompute_corMat <- cbind(cor_with_outcome, expr_corMat)
    precompute_corMat <- rbind(c(1,cor_with_outcome), precompute_corMat)
    cat("Pearson correlation matrix calculation is completed!\n")
    library(parallel)
    #********************************************************
    # Apply PC-simple algorithm around the response variable on NFBC
    #**************************************************************
    pcSimple.fit <- PCSelect_Parallel(y=dfnew[,1],
                                        dm=dfnew[,-1],
                                        method = c("parallel"),
                                        mem.efficient = FALSE,
                                        num_workers = num_workers,
                                        alpha=alpha,
                                        corMat = precompute_corMat,
                                        max_ord=3,
                                        corMethod = "standard",
                                        verbose = TRUE, directed = TRUE)
    # The following codes are used to summarize causally relevant genes
    assoc <- pcSimple.fit$G
    genes <- names(pcSimple.fit$G)
    zMin <- pcSimple.fit$zMin
    res <- cbind(genes,assoc,zMin)
    res_causal <- res[which(assoc==TRUE),]
    colnames(res_causal) <- c("Gene","Assoc","zMin")
    return(res_causal)
}

#' Get cell type specific trait-associated genes for the target traits
#'
#' The following code are used to identify cell type specific associated genes for the target trait
#'
#' @param X_assoc cell type specific gene expression profiles for all genes
#' @param Y the target trait
#' @param PCs_active top PCs used to adjust for population stratificataion
#'
#' @return a result matrix with 5 columns, i.e., Genes, Estimates, Pvalues, Pvalues_BH, Pvalues_Bonferroni
#'
#' @importFrom stats lm
#' 
#' @export
cTWAS <- function(X_assoc,Y,PCs_active){

    #******************************************************************************************************
    # The following code are used to identify cell type specific associated genes for the target trait
    # @funct cTWAS: get cell type specific trait-associated genes for the target traits
    # @param X_assoc: cell type specific gene expression profiles for all genes
    # @param Y: the target trait
    # @param PCs_active: top PCs used to adjust for population stratificataion
    # @return result_mat: a result matrix with 5 columns, i.e., Genes, Estimates, Pvalues, Pvalues_BH, Pvalues_Bonferroni
    # Genes: lists the gene symbols of the causal genes
    # Estimates: The estimated effect size
    # Pvalues: The estimated pvalue
    # Pvalues_BH: The estimated pvalue after correcting for multiple testing by BH
    # Pvalues_Bonferroni: The estimated pvalue after correcting for multiple testing by Bonferroni
    #******************************************************************************************************

    no_genes <- ncol(X_assoc)
    res_mat <- matrix(nrow = no_genes, ncol= 3)
    for (i in 1:no_genes) {
        # cat("This is the adjustment for the",i,"gene",colnames(X_assoc)[i],"\n")
        X <- as.vector(unlist(subset(X_assoc,select=i)))
        input <- cbind(Y,X,PCs_active)
        lm_obj <- lm(Y ~ ., data=input) # use linear regression model, which is much faster than logistic regression model
        #attach(input)
        # lm_obj <- glm( data$outcome ~  gene_expr + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
        res_mat[i,1] <- colnames(X_assoc)[i] # get gene symbol for the target gene
        res_mat[i,2] <- summary(lm_obj)$coeff[2,1] # the estimated effect size
        res_mat[i,3] <- summary(lm_obj)$coeff[2,4] # the estimated Pvalue
    }
    Pvalues_BH <- p.adjust(as.numeric(res_mat[,3]),method ="BH",n = length(res_mat[,3]))
    Pvalue_Bonferroni <- p.adjust(as.numeric(res_mat[,3]),method ="bonferroni",n = length(res_mat[,3]))
    res_mat <- cbind(res_mat,Pvalues_BH,Pvalue_Bonferroni)
    colnames(res_mat) <- c("Genes","Estimates","Pvalues","Pvalues_BH","Pvalues_Bonferroni")
    return(res_mat)
}

#' Detect coefficient change for target genes after adjusting for estimated cell proportions
#' The following code are used to detect coefficient change for target genes after adjusting for estimated cell proportions
#'
#' @param X_assoc tissue-specific gene expression profile for all target genes
#' @param Y the target trait
#' @param PCs_active top PCs used to adjust for population stratificataion
#' @param Cell_prop the estimated cell proportions for the target tissue
#'
#' @return a result matrix with 15 columns, i.e., Genes, Estimates, Pvalues, Pvalues_BH, Pvalues_Bonferroni etc.
#'
#' @importFrom stats lm
#' 
#' @export
coef_change_detection <- function(X_assoc,Y,PCs_active,Cell_prop){

    #**********************************************************************************************************************************
    # The following code are used to detect coefficient change for target genes after adjusting for estimated cell proportions
    # @funct coef_change_detection: detect coefficient change for target genes after adjusting for estimated cell proportions
    # @param X_assoc: tissue-specific gene expression profile for all target genes
    # @param Y: the target trait
    # @param PCs_active: top PCs used to adjust for population stratificataion
    # @param Cell_prop: the estimated cell proportions for the target tissue
    # @return result_mat: a result matrix with 15 columns, i.e.,
    # Genes: lists the gene symbols of the causal genes
    # Estimates(standard TWAS): estimated effect size from standard TWAS analysis
    # Pvalues(standard TWAS): estimated pvalue from standard TWAS analysis
    # SD(standard TWAS): estimated standard deviation from standard TWAS analysis
    # Residual variance(standard TWAS): estimated residual variance from standard TWAS
    # Estimates(conditional TWAS): estimated effect size from conditional TWAS analysis
    # Pvalues(conditional TWAS): estimated pvalue from conditional TWAS analysis
    # SD(conditional TWAS): estimated standard deviation from conditional TWAS analysis
    # Residual variance(conditional TWAS): estimated residual variance from conditional TWAS
    # Coefficient change: coefficient change
    # Pvalue of coefficient change: the estimated pvalue of coefficient change
    # Pvalue Bonferroni(standard TWAS): The estimated pvalue after correcting for multiple testing by Bonferroni from standard TWAS
    # Pvalue BH(standard TWAS): The estimated pvalue after correcting for multiple testing by Bonferroni from standard TWAS
    # Pvalue Bonferroni(conditional TWAS): The estimated pvalue after correcting for multiple testing by Bonferroni from conditional TWAS
    # Pvalue BH(conditional TWAS): The estimated pvalue after correcting for multiple testing by Bonferroni from conditional TWAS
    #***********************************************************************************************************************************

    no_genes <- ncol(X_assoc)
    res_mat <- matrix(nrow = no_genes, ncol= 11 ) # 11 columns 
    for (i in 1:no_genes) {
        X <- as.vector(unlist(subset(X_assoc,select=i)))
        input <- cbind(Y,X,PCs_active)
        lm_obj <- lm(Y ~ ., data=input)
        res_mat[i,1] <- colnames(X_assoc)[i] # get gene symbol for the target gene
        res_mat[i,2] <- summary(lm_obj)$coeff[2,1] # coefficient
        res_mat[i,3] <- summary(lm_obj)$coeff[2,4] # p-value
        res_mat[i,4] <- summary(lm_obj)$coeff[2,2] # standard error
        error_var <- var(lm_obj$resid) # error variance
        res_mat[i,5] <- error_var
        input_cellprop <- cbind(Y,X,PCs_active,Cell_prop)
        lm_obj_cellprop <- lm( Y ~ ., data=input_cellprop )
        res_mat[i,6] <- summary(lm_obj_cellprop)$coeff[2,1] # coefficient
        res_mat[i,7] <- summary(lm_obj_cellprop)$coeff[2,4] # p-value
        res_mat[i,8] <- summary(lm_obj_cellprop)$coeff[2,2] # standard error
        error_var_cellprop <- var(lm_obj_cellprop$resid) # error variance
        res_mat[i,9] <- error_var_cellprop
        coeff_diff <- summary(lm_obj)$coeff[2,1] - summary(lm_obj_cellprop)$coeff[2,1] # coefficient change small model -full model
        res_mat[i,10] <- coeff_diff
        coeff_diff_var <- summary(lm_obj_cellprop)$coeff[2,2]^2 - summary(lm_obj)$coeff[2,2]^2*error_var_cellprop/error_var # formula 15 in Statistical Methods for Comparing Regression Coefficients between Models
        coeff_diff_se <- sqrt(coeff_diff_var)
        coeff_diff_Z <- coeff_diff/coeff_diff_se
        coeff_diff_P <- 2*pnorm(-abs(coeff_diff_Z))
        res_mat[i,11] <- coeff_diff_P # Pvalue for coefficient change test
    }

    Pvalue <- res_mat[,3]
    Pvalue_Bonferroni <- p.adjust(Pvalue,method ="bonferroni",n = length(Pvalue))
    Pvalue_BH <- p.adjust(Pvalue,method ="BH",n = length(Pvalue))

    Pvalue_cellprop <- res_mat[,7]
    Pvalue_cellprop_Bonferroni <- p.adjust(Pvalue_cellprop,method ="bonferroni",n = length(Pvalue_cellprop))
    Pvalue_cellprop_BH <- p.adjust(Pvalue_cellprop,method ="BH",n = length(Pvalue_cellprop))

    res_new <- cbind(res_mat,Pvalue_Bonferroni,Pvalue_BH,Pvalue_cellprop_Bonferroni,Pvalue_cellprop_BH)
    colnames(res_new) <- c("Genes","Estimates(standard TWAS)","Pvalues(standard TWAS)","SD(standard TWAS)","Residual variance(standard TWAS)","Estimates(conditional TWAS)","Pvalues(conditional TWAS)","SD(conditional TWAS)","Residual variance(conditional TWAS)",
                           "Coefficient change","Pvalue of coefficient change","Pvalue Bonferroni(standard TWAS)","Pvalue BH(standard TWAS)","Pvalue Bonferroni(conditional TWAS)","Pvalue BH(conditional TWAS)")
    res_mat <- res_new
    return(res_mat)
}
#' Revise the original gcm.test from R package "GeneralisedCovarianceMeasure"
#' 
#' This code is to revise the original gcm.test from R package "GeneralisedCovarianceMeasure"
#' 
#' @param X A(nxp)-dimensional matrix (or data frame) with n observations of p variables.
#' @param Y A(nxp)-dimensional matrix (or data frame) with n observations of p variables.
#' @param Z A(nxp)-dimensional matrix (or data frame) with n observations of p variables.
#' @param alpha Significance level of the test.
#' @param regr.method A string indicating the regression method that is used. Currently implemented are "gam", "xgboost", "kernel.ridge". The regression is performed only if not both resid.XonZ and resid.YonZ are set to NULL.
#' @param regr.pars Some regression methods require a list of additional options
#' @param plot.residuals A Boolean indicating whether some plots should be shown.
#' @param nsim An integer indicating the number of bootstrap samples used to approximate the null distribution of the test statistic.
#' @param resid.XonZ It is possible to directly provide the residuals instead of performing a regression. If set to NULL, the regression method specified in regr.method is used.
#' @param resid.YonZ It is possible to directly provide the residuals instead of performing a regression. If set to NULL, the regression method specified in regr.method is used.
#' 
#' @return a list containing P-value, R, test.statistic and whether reject
#' 
#' @import GeneralisedCovarianceMeasure
gcm_test <- function(X = NULL, Y = NULL, Z = NULL, alpha = 0.05, regr.method = "xgboost", 
                     regr.pars = list(), plot.residuals = FALSE, nsim = 499L, 
                     resid.XonZ = NULL, resid.YonZ = NULL) 
{

    # Written by yinly, Nov 2, 2022
    # Description: revise the original gcm.test from R package "GeneralisedCovarianceMeasure"
    # output calculated product of residual for each subject, which later on could be used to detect the signficance
    # of statistic change between two gcm tests

    if (is.null(Z)) {
        if (is.null(resid.XonZ)) {
            resid.XonZ <- X
        }
        if (is.null(resid.YonZ)) {
            resid.YonZ <- Y
        }
    }
    else {
        if (is.null(resid.XonZ)) {
            if (is.null(X)) {
                stop("Either X or resid.XonZ must be provided.")
            }
            resid_func <- function(V) comp.resids(V, Z, regr.pars, 
                regr.method)
            if (is.matrix(X)) {
                resid.XonZ <- apply(X, 2, resid_func)
            }
            else {
                resid.XonZ <- resid_func(X)
            }
        }
        if (is.null(resid.YonZ)) {
            if (is.null(Y)) {
                stop("Either Y or resid.YonZ must be provided.")
            }
            resid_func <- function(V) comp.resids(V, Z, regr.pars, 
                regr.method)
            if (is.matrix(Y)) {
                resid.YonZ <- apply(Y, 2, resid_func)
            }
            else {
                resid.YonZ <- resid_func(Y)
            }
        }
    }
    nn <- NA
    if (NCOL(resid.XonZ) > 1 || NCOL(resid.YonZ) > 1) {
        d_X <- NCOL(resid.XonZ)
        d_Y <- NCOL(resid.YonZ)
        nn <- NROW(resid.XonZ)
        R_mat <- rep(resid.XonZ, times = d_Y) * as.numeric(as.matrix(resid.YonZ)[, 
            rep(seq_len(d_Y), each = d_X)])
        dim(R_mat) <- c(nn, d_X * d_Y)
        R_mat <- t(R_mat)
        R_mat <- R_mat/sqrt((rowMeans(R_mat^2) - rowMeans(R_mat)^2))
        test.statistic <- max(abs(rowMeans(R_mat))) * sqrt(nn)
        test.statistic.sim <- apply(abs(R_mat %*% matrix(rnorm(nn * 
            nsim), nn, nsim)), 2, max)/sqrt(nn)
        p.value <- (sum(test.statistic.sim >= test.statistic) + 
            1)/(nsim + 1)
        if (plot.residuals) {
            par(mfrow = c(NCOL(resid.XonZ), NCOL(resid.YonZ)))
            for (ii in 1:NCOL(resid.XonZ)) {
                for (jj in 1:NCOL(resid.YonZ)) {
                  plot(resid.XonZ[, ii], resid.YonZ[, jj], main = "scatter plot of residuals")
                }
            }
            par(mfrow = c(1, 1))
        }
    }
    else {
        nn <- ifelse(is.null(dim(resid.XonZ)), length(resid.XonZ), dim(resid.XonZ)[1])
        R <- resid.XonZ * resid.YonZ
        R.sq <- R^2
        meanR <- mean(R)
        test.statistic <- sqrt(nn) * meanR/sqrt(mean(R.sq) - 
            meanR^2)
        p.value <- 2 * pnorm(abs(test.statistic), lower.tail = FALSE)
        if (plot.residuals) {
            plot(resid.XonZ, resid.YonZ, main = "scatter plot of residuals")
        }
    }
    return(list(p.value = p.value, R = R, test.statistic = test.statistic, 
           reject = (p.value < alpha)))
}

#' Calculate the covariance between the calculated GCMs from two different GCM tests
#' 
#' The following code is used to calculate the covariance between the calculated GCMs from two different GCM tests.
#' 
#' @param R1 derived residual product for the 1st test
#' @param R2 derived residual product for the 2nd test
#' @param N sample size
#' 
#' @return calculated covariance
stat_change_detection <- function(R1, R2, N){

    # Description: quantify the significance of statistics change derived from GCM
    # R1: derived residual product for the 1st test
    # R2: derived residual product for the 2nd test
    # N: sample size

    se_R1 <- sqrt(var(R1))/sqrt(N)
    se_R2 <- sqrt(var(R2))/sqrt(N)
    tn_cov <- 1/(se_R1*se_R2) * cov(R1,R2)/N
    return(tn_cov)
}

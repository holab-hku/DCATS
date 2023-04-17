
#' A beta-binomial regression test with similarity based bootstrapping
#'
#' A GLM test with binomial distribution. In order to estimate the variance of
#' the weight, a boostrapping based on the composition similarity is performed.
#'
#' @param counts1 A matrix of compsition sizes (n_rep1, n_cluster) for each
#'   replicate in each cluster for codition 1 as case
#' @param counts2 A matrix of compsition sizes (n_rep2, n_cluster) for each
#'   replicate in each cluster for codition 2 as control
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#'   similarity matrix between cluster group pair. The order of cluster should
#'   be consistent with those in `counts1` and `counts2`
#' @param pseudo_count A pseudo count to add for counts in all cell types.
#'   Default NULL means 0 except if a cell type is emplty in one condition,
#'   otherwise pseudo_count will be: 0.01 * rowMeans for each condition
#' @param n_samples An integer for number samples in sampling for estimating the
#'   variance of the weights
#'
#' @return a vector of significance p values for each cluster
#'
#' @importFrom matrixStats colSds
#' 
#'
#' @examples
#' K <- 2
#' totals1 = c(100, 800, 1300, 600)
#' totals2 = c(250, 700, 1100)
#' diri_s1 = rep(1, K) * 20
#' diri_s2 = rep(1, K) * 20
#' simil_mat = create_simMat(K, confuse_rate=0.2)
#' sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#' #dcats_betabin(sim_dat[[1]], sim_dat[[2]], simil_mat, n_samples = 100)
#' 

dcats_betabin <- function(counts1, counts2, similarity_mat=NULL, n_samples=50,
                          pseudo_count=NULL) {
    ## Check counts1 and counts2 shape
    if (length(counts1) == 1 || is.null(dim(counts1)) ||
        length(dim(counts1)) < 2) {
        counts1 = matrix(counts1, nrow=1)
    }
    if (length(counts2) == 1 || is.null(dim(counts2)) ||
        length(dim(counts2)) < 2) {
        counts2 = matrix(counts2, nrow=1)
    }

    prop1 <- counts1 / rowSums(counts1)
    prop2 <- counts2 / rowSums(counts2)

    ## using estimated the latent cell counts
    counts1_latent = counts1
    counts2_latent = counts2
    if(!is.null(similarity_mat)) {
        for (i in seq_len(nrow(counts1))) {
            counts1_latent[i, ] <- sum(counts1[i, ]) *
                multinom_EM(counts1[i, ], similarity_mat, verbose = FALSE)$mu
        }
        for (i in seq_len(nrow(counts2))) {
            counts2_latent[i, ] <- sum(counts2[i, ]) *
                multinom_EM(counts2[i, ], similarity_mat, verbose = FALSE)$mu
        }
    }
    
    ## number of cell types
    K <- ncol(counts1)
    if (is.null(similarity_mat)) {
        n_samples <- 1
    }

    if (!is.null(n_samples) && !is.null(similarity_mat)) {
        counts1_use <- matrix(0, nrow(counts1) * n_samples, K)
        counts2_use <- matrix(0, nrow(counts2) * n_samples, K)
        for (i in seq_len(nrow(counts1))) {
            idx <- seq((i - 1) * n_samples + 1, i * n_samples)
            for (j in seq_len(K)) {
                counts1_use[idx, ] <- (
                    counts1_use[idx, ] + t(rmultinom(n_samples,
                                                     counts1_latent[i, j],
                                                     similarity_mat[j, ])))
            }
        }
        for (i in seq_len(nrow(counts2))) {
            idx <- seq((i - 1) * n_samples + 1, i * n_samples)
            for (j in seq_len(K)) {
                counts2_use[idx, ] <- (
                    counts2_use[idx, ] + t(rmultinom(n_samples,
                                                     counts2_latent[i, j],
                                                     similarity_mat[j, ])))
            }
        }
    } else{
        counts1_use <- counts1
        counts2_use <- counts2
    }

    # adding pseudo counts
    if (is.null(pseudo_count)) {
        if (any(colMeans(counts1) == 0) || any(colMeans(counts2) == 0) ) {
            warning(paste("Empty cell type exists in at least one conidtion;",
                        "adding replicate & condition specific pseudo count:"))
            counts1_use <- counts1_use + 1
            counts2_use <- counts2_use + 1
        }
    } else {
        counts1_use = counts1_use + pseudo_count
        counts2_use = counts2_use + pseudo_count
    }

    ## Binomial regression for each sampling
    LR_val <- matrix(NA, n_samples, K)
    coeffs_val <- matrix(NA, n_samples, K)
    coeffs_err <- matrix(NA, n_samples, K)
    intercept_val <- matrix(NA, n_samples, K)
    intercept_err <- matrix(NA, n_samples, K)
    total_all <- c(rowSums(counts1_use), rowSums(counts2_use))
    label_all <- c(rep(1, nrow(counts1_use)), rep(0, nrow(counts2_use)))
    for (ir in seq_len(n_samples)) {           ## for each sampling
        idx <- seq(1, length(total_all), n_samples) + ir - 1
        for (i in seq_len(K)) {                ## for each cluster
            n1 <- c(counts1_use[, i], counts2_use[, i])[idx]
            df <- data.frame(n1 = n1, n2 = total_all[idx] - n1,
                             label = label_all[idx])
            ## betabin GLM
            fm1 <- aod::betabin(cbind(n1, n2) ~ label, ~ 1, data = df)
            if (any(is.na(fm1@varparam))) {
                message("Existing NA in fm1.")
            } else {
                coeffs_err[ir, i] <-fm1@varparam[2, 2]
                }
            coeffs_val[ir, i] <- fm1@param[2]
            intercept_val[ir, i] <- fm1@param[1]
            intercept_err[ir, i] <- fm1@varparam[1, 1]
            
            fm0 <- aod::betabin(cbind(n1, n2) ~ 1, ~ 1, data = df)
            LR_val[ir, i] <- fm0@dev - fm1@dev
            }
    }

    ## Averaging the coeficients errors
    if (is.null(n_samples) || is.null(similarity_mat) || n_samples == 1) {
        coeff_val_mean <- colMeans(coeffs_val)
        coeff_err_pool <- colMeans(coeffs_err**2)
        intercept_val_mean <- colMeans(intercept_val)
        intercept_err_pool <- colMeans(intercept_err**2)
    } else{
        coeff_val_mean <- colMeans(coeffs_val)
        coeff_err_pool <- colMeans(coeffs_err**2) +
            matrixStats::colSds(coeffs_val) +
            matrixStats::colSds(coeffs_val) / n_samples

        intercept_val_mean <- colMeans(intercept_val)
        intercept_err_pool <- colMeans(intercept_err**2) +
            matrixStats::colSds(intercept_val) +
            matrixStats::colSds(intercept_val) / n_samples
    }

    # p values with Ward test: https://en.wikipedia.org/wiki/Wald_test
    pvals <- pnorm(-abs(coeff_val_mean) / sqrt(coeff_err_pool))  * 2

    LR_median = robustbase::colMedians(LR_val)
    LRT_pvals <-  pchisq(LR_median, df=1, lower.tail = FALSE, log.p = FALSE)

    data.frame(
        "prop1_mean" = colMeans(prop1),
        "prop1_std"  = matrixStats::colSds(prop1),
        "prop2_mean" = colMeans(prop2),
        "prop2_std"  = matrixStats::colSds(prop2),
        "coeff_mean" = coeff_val_mean,
        "coeff_std"  = sqrt(coeff_err_pool),
        "intecept_mean" = intercept_val_mean,
        "intecept_std"  = sqrt(intercept_err_pool),
        "pvals" = pvals,
        "LRT_pvals" = LRT_pvals,
        row.names = colnames(counts1)
    )
}



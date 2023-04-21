#' A Generalised linear model based likelihood ratio testing
#'
#' GLM supports both beta-binomial and negative binomial from aod package.
#'
#' @param count_mat A matrix of composition sizes (n_sample, n_cluster) for each
#'   cluster in each sample
#' @param design_mat A matrix or a data frame of testing candidate factors 
#'   (n_sample, n_factor) with same sample order as count_mat. All factors 
#'   should be continous and categorical with only two levels. 
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#'   similarity matrix between cluster group pair. The order of cluster should
#'   be consistent with those in `count_mat`.
#' @param pseudo_count A pseudo count to add for counts in all cell types
#'   Default NULL means 0 except if a cell type is empty in one condition,
#'   otherwise pseudo_count will be: 0.01 * rowMeans for each condition
#' @param base_model A string value: `NULL` for 1 factor vs NULL factor testing;
#'   `FULL` for FULL factors vs n-1 factors testing.
#' @param fix_phi A numeric used to provided a fixed phi value for the GLM 
#'   for all cell types
#' @param reference A vector of characters indicating which cell types are used 
#' as reference for normalization. `NULL` indicates using total count for 
#' normalization.
#'
#' @return a list of significance p values for each cluster
#'
#' @export
#' @importFrom matrixStats colSds
#' @import stats
#'
#' @examples
#' K <- 3
#' totals1 = c(100, 800, 1300, 600)
#' totals2 = c(250, 700, 1100)
#' diri_s1 = rep(1, K) * 20
#' diri_s2 = rep(1, K) * 20
#' simil_mat = DCATS::create_simMat(K, confuse_rate=0.2)
#' sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#' sim_count = rbind(sim_dat$numb_cond1, sim_dat$numb_cond2)
#' sim_design = data.frame(condition = c("g1", "g1", "g1", "g1", "g2", "g2", "g2"), 
#' gender = sample(c("Female", "Male"), 7, replace = TRUE))
#' ## Using 1 factor vs NULL factor testing
#' dcats_GLM(sim_count, sim_design, similarity_mat = simil_mat)
#' ## Using full factors vs n-1 factors testing with intercept term
#' dcats_GLM(sim_count, sim_design, similarity_mat = simil_mat, base_model='FULL')
#' ## Fix phi
#' dcats_GLM(sim_count, sim_design, similarity_mat = simil_mat, fix_phi = 1/61)
#' ## Specify reference cell type
#' colnames(sim_count) <- c("celltypeA", "celltypeB", "celltypeC")
#' 

dcats_GLM <- function (count_mat, design_mat, similarity_mat = NULL, pseudo_count = NULL, 
                       base_model = "NULL", fix_phi = NULL, reference = NULL) 
{
  coeffs <- matrix(NA, ncol(count_mat), ncol(design_mat))
  coeffs_err <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LR_vals <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_pvals <- matrix(NA, ncol(count_mat), ncol(design_mat))
  pvals <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_fdr <- matrix(NA, ncol(count_mat), ncol(design_mat))
  if (is.null(colnames(count_mat))) 
    colnames(count_mat) <- paste0("cell_type_", seq(ncol(count_mat)))
  if (is.null(colnames(design_mat))) 
    colnames(design_mat) <- paste0("factor_", seq(ncol(design_mat)))
  rownames(LR_vals) <- rownames(LRT_pvals) <- rownames(pvals) <- rownames(LRT_fdr) <- rownames(coeffs) <- rownames(coeffs_err) <- colnames(count_mat)
  colnames(LR_vals) <- colnames(LRT_pvals) <- colnames(pvals) <- colnames(LRT_fdr) <- colnames(coeffs) <- colnames(coeffs_err) <- colnames(design_mat)
  count_use <- count_mat
  if (!is.null(similarity_mat)) {
    for (i in seq_len(nrow(count_mat))) {
      count_use[i, ] <- sum(count_mat[i, ]) * multinom_EM(count_mat[i, 
      ], similarity_mat)$mu
    }
  }
  K <- ncol(count_mat)
  if (is.null(pseudo_count)) {
    if (any(colMeans(count_mat) == 0)) {
      warning(paste("Empty cell type exists in at least one conidtion;", 
                  "adding replicate & condition specific pseudo count:"))
      count_use <- count_use + 1
    }
  }
  else {
    count_use <- count_use + pseudo_count
  }
  count_use <- round(count_use)
  n_samples <- 1
  for (k in seq_len(ncol(design_mat))) {
    sub_LR_val <- matrix(NA, n_samples, K)
    sub_coeffs_val <- matrix(NA, n_samples, K)
    sub_coeffs_err <- matrix(NA, n_samples, K)
    for (ir in seq_len(n_samples)) {
      idx <- seq(1, nrow(count_use), n_samples) + ir - 1
      for (m in seq_len(ncol(count_use))) {
        if (is.null(reference)) {
          df_use <- data.frame(n1 = count_use[, m], total = rowSums(count_use))[idx, 
          ]
          df_use$ref_count <- df_use$total - df_use$n1
        }
        else {
          if (length(reference) == 1) {df_use <- data.frame(n1 = count_use[, m], ref_count = count_use[, 
                                                                                                       reference])[idx, ]} else {df_use <- data.frame(n1 = count_use[, m], ref_count = rowSums(count_use[, 
                                                                                                                                                                                                         reference]))[idx, ]}
        }
        df_use <- cbind(df_use, design_mat)
        df_tmp <- df_use[!is.na(design_mat[, k]), ]
        if (base_model == "NULL" | ncol(design_mat) == 
            1) {
          formula_fm0 <- as.formula("cbind(n1, ref_count) ~ 1")
          formula_fm1 <- as.formula(paste0("cbind(n1, ref_count)", 
                                           "~ 1+", colnames(design_mat)[k], sep = ""))
        }
        else if (base_model == "FULL") {
          fm0_right <- paste(colnames(design_mat)[-k], 
                             collapse = " + ")
          fm1_right <- paste(colnames(design_mat), collapse = " + ")
          formula_fm0 <- as.formula(paste0("cbind(n1, ref_count)", 
                                           " ~ 1 + ", fm0_right, sep = ""))
          formula_fm1 <- as.formula(paste0("cbind(n1, ref_count)", 
                                           " ~ 1 + ", fm1_right, sep = ""))
        }
        fm0 <- aod::betabin(formula_fm0, ~1, data = df_tmp, 
                            warnings = FALSE)
        fm1 <- aod::betabin(formula_fm1, ~1, data = df_tmp, 
                            warnings = FALSE)
        if (!is.null(fix_phi)) {
          fm0 <- aod::betabin(formula_fm0, ~1, data = df_tmp, 
                              warnings = FALSE, fixpar = list(fm0@nbpar, 
                                                              fix_phi))
          fm1 <- aod::betabin(formula_fm1, ~1, data = df_tmp, 
                              warnings = FALSE, fixpar = list(fm1@nbpar, 
                                                              fix_phi))
        }
        if (length(fm1@varparam) < 4 || is.na(fm1@varparam[2, 
                                                           2])) {
          next
        }
        sub_LR_val[ir, m] <- fm0@dev - fm1@dev
        parID <- grep(colnames(design_mat)[k], names(fm1@param))
        if (length(parID) > 1) 
          stop("Please check the design matrix, make sure all factors are continous or categorical with only two levels.")
        sub_coeffs_val[ir, m] <- fm1@param[parID]
        if (is.null(fix_phi)) 
          sub_coeffs_err[ir, m] <- fm1@varparam[parID, 
                                                parID]
      }
    }
    coeff_val_mean <- colMeans(sub_coeffs_val, na.rm = TRUE)
    if (is.null(fix_phi)) {
      if (is.null(n_samples) || is.null(similarity_mat) || 
          n_samples == 1) {
        sub_coeff_err_pool <- colMeans(sub_coeffs_err^2, 
                                       na.rm = TRUE)
      }
      else {
        sub_coeff_err_pool <- colMeans(sub_coeffs_err^2, 
                                       na.rm = TRUE) + matrixStats::colSds(sub_coeffs_val) + 
          matrixStats::colSds(sub_coeffs_val)/n_samples
      }
      pvals[, k] <- pnorm(-abs(coeff_val_mean)/sqrt(sub_coeff_err_pool)) * 
        2
      coeffs_err[, k] <- sqrt(sub_coeff_err_pool)
    }
    LR_median <- robustbase::colMedians(sub_LR_val, na.rm = TRUE)
    LR_vals[, k] <- LR_median
    LRT_pvals[, k] <- pchisq(LR_median, df = 1, lower.tail = FALSE, 
                             log.p = FALSE)
    coeffs[, k] <- coeff_val_mean
  }
  LRT_fdr[, ] <- p.adjust(LRT_pvals, method = "fdr")
  res <- list(ceoffs = coeffs, coeffs_err = coeffs_err, LR_vals = LR_vals, 
              LRT_pvals = LRT_pvals, fdr = LRT_fdr)
  res
}

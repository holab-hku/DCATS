#' Calculate a global phi for all cell types
#'
#' Assuming all cell types share the same phi. This global phi can be calculate 
#' by pooling all cell types together to fit a beta binomial distribution.
#'
#' @param count_mat A matrix of composition sizes (n_sample, n_cluster) for each
#'   cluster in each sample
#' @param design_mat A matrix of testing candidate factors (n_sample, n_factor)
#'   with same sample order as count_mat
#'   
#' @return A number indicating a global phi for all cell types
#'   
#' @export
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
#' sim_design = data.frame(condition = c("g1", "g1", "g1", "g1", "g2", "g2", "g2"))
#' phi = DCATS::getPhi(sim_count, sim_design)


getPhi <- function(count_mat, design_mat){
  K <- ncol(count_mat)
  S <- nrow(count_mat)
  total <- rep(rowSums(count_mat), K)
  n1 <- matrix(count_mat, ncol = 1)
  design_use <- do.call("rbind", replicate(K, design_mat, simplify = FALSE))
  celltype_idx <- matrix(0, S*K, K)
  for (type in seq(1,K)) {
    celltype_idx[(1:S)+S*(type-1),type] = 1
  }
  df_use <- data.frame(total, n1, design_use, celltype_idx[,-1])
  fm <- aod::betabin(cbind(n1, total-n1) ~ ., ~ 1, data = df_use, warnings = FALSE) # make intercept as 0
  return(fm@random.param)
}

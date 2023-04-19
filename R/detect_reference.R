#' Calculate a global phi for all cell types
#'
#' Assuming all cell types share the same phi. This global phi can be calculate 
#' by pooling all cell types together to fit a beta binomial distribution.
#'
#' @param count_mat A matrix of composition sizes (n_sample, n_cluster) for each
#'   cluster in each sample.
#' @param design_mat A matrix or a data frame of testing candidate factors 
#'   (n_sample, n_factor) with same sample order as count_mat. All factors 
#'   should be continous and categorical with only two levels. 
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#'   similarity matrix between cluster group pair. The order of cluster should
#'   be consistent with those in `count_mat`.
#' @param fix_phi A numeric used to provided a fixed phi value for the GLM 
#'   for all cell types.
#'   
#' @return A data frame with ordered cell types and their p-value. Cell types 
#' are ordered by their p-values. The order indicating how they are recommended
#' to be selected as reference cell types.
#' 
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
#' sim_design = data.frame(condition = c("g1", "g1", "g1", "g1", "g2", "g2", "g2"), 
#' gender = sample(c("Female", "Male"), 7, replace = TRUE))
#' ## Using 1 factor vs NULL factor testing
#' detect_reference(sim_count, sim_design)

detect_reference <- function(count_mat, design_mat, similarity_mat = NULL, fix_phi = NULL){
  res <- dcats_GLM(count_mat, design_mat, similarity_mat, fix_phi = fix_phi)
  resDF <- data.frame(celltype = rownames(res$LRT_pval), pval = as.vector(apply(res$LRT_pvals, 1, min)))
  resDF <- resDF[resDF$pval > 0.1,]
  resDF <- resDF[order(-resDF$pval),]
  if (nrow(resDF) < 2)
    return("No suitable reference cell type detected.")
  ## calculate reference group proportion
  if (is.null(colnames(count_mat))) 
    colnames(count_mat) <- paste0("cell_type_", seq(ncol(count_mat)))
  typeSum <- colSums(count_mat)
  for (celltype_num in seq(2,nrow(resDF))) {
    if (sum(typeSum[resDF$celltype[1:celltype_num]])/sum(typeSum) > 0.25) {
      message(paste0("Please use at least ", as.character(celltype_num), " cell types as the reference group."))
      return(resDF$celltype)
    }
  }
}

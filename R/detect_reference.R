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

detect_reference = function(count_mat, design_mat, similarity_mat = NULL, fix_phi = NULL){
  res = dcats_GLM(count_mat, design_mat, similarity_mat, fix_phi = fix_phi)
  resDF = data.frame(celltype = rownames(res$LRT_pval), pval = as.vector(res$LRT_pvals))
  resDF = resDF[order(-resDF$pval),]
  return(resDF)
}

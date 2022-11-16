#' Generate similarity matrix with uniform confusion rate to none-self clusters
#'
#' Create a similarity matrix assuming the misclassification error distribute
#' uniformly in all clusters
#'
#' @param K A integer for number of cluster
#' @param confuse_rate A float for confusion rate, uniformly to none-self
#'   clusters
#'
#' @return a similarity matrix with uniform confusion with other cluster
#'
#' @export
#'
#' @examples
#' create_simMat(4, 0.1)
#' 

create_simMat <- function(K, confuse_rate) {
  diag(K) * (1 - confuse_rate) + confuse_rate * (1 - diag(K)) / (K - 1)
}
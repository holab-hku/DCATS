
#' Composition size simulator
#'
#' Directly simulating the composition size from a Dirichlet-Multinomial
#' distribution with replicates for two conditions.
#'
#' @param totals_cond1 A vector of integers (n_rep1, ) for the total samples in
#'   each replicate in codition 1
#' @param totals_cond2 A vector of integers (n_rep2, ) for the total samples in
#'   each replicate in codition 2
#' @param dirichlet_s1 A vector of floats (n_cluster, ) for the composition
#'   concentration in condition 1
#' @param dirichlet_s2 A vector of floats (n_cluster, ) for the composition
#'   concentration in condition 2
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#'   similarity matrix between each cluster pair
#'
#' @return a list of two matrices for composition sizes in each replicate and
#'   each cluster in both conditions.
#'
#' @export
#' @importFrom MCMCpack rdirichlet
#'
#'
#' @examples
#' K <- 2
#' totals1 = c(100, 800, 1300, 600)
#' totals2 = c(250, 700, 1100)
#' diri_s1 = rep(1, K) * 20
#' diri_s2 = rep(1, K) * 20
#' confuse_rate = 0.2
#' simil_mat = create_simMat(2, 0.2)
#' sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#' 
simulator_base <- function(totals_cond1, totals_cond2,
                           dirichlet_s1, dirichlet_s2,
                           similarity_mat=NULL) {
    K <- length(dirichlet_s1)       # number of clusters
    n_rep1 <- length(totals_cond1)  # number of replicates in condition 1
    n_rep2 <- length(totals_cond2)  # number of replicates in condition 2

    prop_cond1 <- matrix(0, n_rep1, K)
    prop_cond2 <- matrix(0, n_rep2, K)
    numb_cond1 <- matrix(0, n_rep1, K)
    numb_cond2 <- matrix(0, n_rep2, K)
    for (i in seq_len(n_rep1)) {
        prop_cond1[i, ] <- MCMCpack::rdirichlet(1, dirichlet_s1)
        numb_cond1[i, ] <- rmultinom(1, totals_cond1[i], prop_cond1[i, ])
    }
    for (i in seq_len(n_rep2)) {
        prop_cond2[i, ] <- MCMCpack::rdirichlet(1, dirichlet_s2)
        numb_cond2[i, ] <- rmultinom(1, totals_cond2[i], prop_cond2[i, ])
    }

    if (is.null(similarity_mat)) {
        similarity_mat <- diag(K)
    }

    ## sampling based on similarity matrix
    numb_out_cond1 <- matrix(0, n_rep1, K)
    numb_out_cond2 <- matrix(0, n_rep2, K)
    for (i in seq_len(n_rep1)) {
        for (j in seq_len(K)) {
            numb_out_cond1[i, ] <- (numb_out_cond1[i, ] +
                                        rmultinom(1, numb_cond1[i, j],
                                                  similarity_mat[j, ])[, 1])
        }
    }

    for (i in seq_len(n_rep2)) {
        for (j in seq_len(K)) {
            numb_out_cond2[i, ] <- (numb_out_cond2[i, ] +
                                        rmultinom(1, numb_cond2[i, j],
                                                  similarity_mat[j, ])[, 1])
        }
    }

    list("numb_cond1" = numb_out_cond1,
         "numb_cond2" = numb_out_cond2,
         "true_cond1" = numb_cond1,
         "true_cond2" = numb_cond2)
}


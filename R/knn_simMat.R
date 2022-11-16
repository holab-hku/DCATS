
#' Calculate stochastic transition matrix between clusters from a KNN connection
#' matrix
#'
#' The transition probability from cluster i to j is the fraction of neighbours
#' of all samples in cluster i that belongs to cluster j. Note, this matrix is
#' asymmetric, so as the input KNN connection matrix.
#'
#' @param KNN_matrix a sparse binary matrix with size (n_sample, n_sample).
#'   x_ij=1 means sample j is a neighbour of sample i. As definition, we expect
#'   sum(KNN_matrix) = n_sample * K, where K is the number neighbours.
#' @param clusters a (n_sample, ) vector of cluster id for each sample.
#' 
#' @export
#'
#' @examples
#' data(simulation)
#' knn_mat = knn_simMat(simulation$knnGraphs, simulation$labels)
#' 
#' @return a similarity matrix calculated based on the knn graph.

knn_simMat <- function(KNN_matrix, clusters){
    clusters_uniq <- unique(clusters)
    n_cluster <- length(clusters_uniq)

    simi_mat <- matrix(0, n_cluster, n_cluster)
    rownames(simi_mat) <- colnames(simi_mat) <- clusters_uniq
    for (i in seq_len(n_cluster)) {
        idx_i <- clusters == clusters_uniq[i]
        for (j in seq_len(n_cluster)) {
            idx_j <- clusters == clusters_uniq[j]
            simi_mat[i, j] <- sum(KNN_matrix[idx_i, ][, idx_j])
        }
    }
    simi_mat / rowSums(simi_mat)
}

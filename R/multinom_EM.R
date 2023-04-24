#' An EM algorithm to fit a multinomial with maximum likelihood
#'
#' @param X A vector of commopent sizes
#' @param simMM A matrix of floats (n_cluster, n_cluster) for the similarity
#'   matrix between clusters. simMM[i,j] means the proportion of cluster i will
#'   be assigned to cluster j, hence colSums(simMM) are ones.
#' @param max_iter integer(1). number of maximum iterations
#' @param min_iter integer(1). number of minimum iterations
#' @param logLik_threshold A float. The threshold of logLikelihood increase for
#'   detecting convergence
#'
#' @return a list containing \code{mu}, a vector for estimated latent proportion
#'   of each cluster, \code{logLik}, a float for the estimated log likelihood,
#'   \code{simMM}, the input of simMM, code{X}, the input of X, \code{X_prop},
#'   the proportion of clusters in the input X, \code{predict_X_prop}, and the
#'   predicted proportion of clusters based on mu and simMM.
#'
#' @export
#'
#' @examples
#' X = c(100, 300, 1500, 500, 1000)
#' simMM = create_simMat(5, confuse_rate=0.2)
#' multinom_EM(X, simMM)
#' 
multinom_EM <- function(X, simMM, min_iter=10, max_iter=1000, logLik_threshold=1e-2) {
    # Be very careful on the shape of simMM; rowSums(simMM) = 1
    K <- ncol(simMM)

    # initialization
    mu <- sample(K)
    mu <- mu / sum(mu)
    Z <- matrix(NA, K, K)
    logLik_old <- logLik_new <- log(mu %*% simMM) %*% X

    for (it in seq_len(max_iter)) {
        ## E step: expectation of count from each component

        for (i in seq(K)) {
            for (j in seq(K)){
                Z[i, j] <- simMM[i, j] * mu[i] / sum(mu * simMM[, j])
            }
        }

        ## M step: maximizing likelihood

        ## v2
        mu <- c(Z %*% X)
        mu <- mu / sum(mu)

        ## Check convergence
        logLik_new <- log(mu %*% simMM) %*% X
        if (it > min_iter && logLik_new - logLik_old < logLik_threshold) {
            break
        } else {
            logLik_old <- logLik_new
        }

    }

    ## return values
    list("mu" = mu, "logLik" = logLik_new, "simMM" = simMM, "X" = X, 
         "X_prop" = X / sum(X), "predict_X_prop" = mu %*% simMM)
}

# matrix(seq(6), 2, 3) * c(1, 2)

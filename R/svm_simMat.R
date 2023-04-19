
#' Calculate stochastic transition matrix between clusters from a data frame 
#' including information about clustering
#'
#' The transition probability from cluster i to j is calculated based on the 
#' information used to cluster cells. It is estimated by the misclassification 
#' rate from cluster i to j comparing the original labels with the labels predicted 
#' by support vector machine with 5-fold cross validation.
#'
#' @param dataframe a data frame contains the information used for clustering 
#' and the original label of each cell. The original labels should have the 
#' column name `clusterRes`.
#' 
#' @export
#' @examples
#' data(Kang2017)
#' svm_mat = svm_simMat(Kang2017$svmDF)
#'
#' @return a similarity matrix estimated by 5-fold cross validation support 
#' vector machine.

svm_simMat <- function(dataframe){
  dataframe$clusterRes <- as.factor(dataframe$clusterRes)
  row_num <- nrow(dataframe)
  pred <- rep(NA, row_num)
  idx1 <- sample(seq(1,row_num), round(row_num/5))
  idx2 <- sample(seq(1,row_num)[-idx1], round(row_num/5))
  idx3 <- sample(seq(1,row_num)[-c(idx1, idx2)], round(row_num/5))
  idx4 <- sample(seq(1,row_num)[-c(idx1, idx2, idx3)], round(row_num/5))
  idx5 <- seq(1,row_num)[-c(idx1, idx2, idx3, idx4)]
  idxV <- list(idx1, idx2, idx3, idx4, idx5)
  for (i in seq(1,5)){
    svmfit <- e1071::svm(clusterRes ~ ., data = dataframe[-idxV[[i]],], kernel = "radial")
    pred[idxV[[i]]] <- as.character(aod::predict(svmfit, subset(dataframe, select = -clusterRes)[idxV[[i]],]))
  }
  conf.mat <- table(dataframe$clusterRes, as.factor(pred))
  simil_mat <- t(t(conf.mat)/apply(conf.mat,2,sum))
  return(simil_mat)
}


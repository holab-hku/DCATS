#' Count matrices of intestinal epithelial scRNA-seq data from three conditions
#'
#' A data containing the count matrices, the similarity matrix and other 
#' variables used to generate the similarity matrix 
#' from intestinal epithelial single cell RNA sequencing data with 
#' three condition
#'
#' @format A list with 7 items:
#' \describe{
#'   \item{count_ctrl}{the count matrix for the control group}
#'   \item{count_Hpoly3}{the count matrix for three days after H.polygyrus infection}
#'   \item{count_Hpoly10}{the count matrix for ten days after H.polygyrus infection}
#'   \item{count_Salma}{the count matrix for two days after Salmonella infection}
#'   \item{svm_mat}{the similarity matrix}
#'   \item{svmDF}{the data frame used to calculate the similarity matrix}
#'   \item{source}{the source of this dataset}
#' }
#' @source \url{https://www.nature.com/articles/nature24489}
#' @examples 
#' library(DCATS)
#' data(Haber2017)
"Haber2017"

#' Count matrices of 8 pooled lupus patient samples within two conditions
#'
#' A data containing the count matrices, the similarity matrix and other 
#' variables used to generate the similarity matrix 
#' from single cell RNA sequencing data of 8 pooled lupus patient samples within 
#' two conditions
#'
#' @format A list with 5 items:
#' \describe{
#'   \item{count_ctrl}{the count matrix for three days after H.polygyrus infection}
#'   \item{count_stim}{the count matrix for ten days after H.polygyrus infection}
#'   \item{svm_mat}{the simularity matrix}
#'   \item{svmDF}{the data frame used to calculate the similarity matrix.}
#'   \item{source}{the source of this dataset}
#' }
#' @source \url{https://www.nature.com/articles/nbt.4042}
"Kang2017"

#' Count matrix and metadata of a large COVID-19 scRNA-seq data cohort
#'
#' A data containing the count matrix, metadata from a large COVID-19 cohort
#'
#' @format A list with 7 items:
#' \describe{
#'   \item{countM}{the count matrix for all samples comming from different condition}
#'   \item{designM}{the corresponding metadata related to each sample}
#'   \item{source}{the source of this dataset}
#' }
#' @source \url{https://www.nature.com/articles/nature24489}
"Ren2021"

#' Simulated dataset with two conditions
#'
#' A data containing the count matrices, the similarity matrix and other 
#' variables used to generate the similarity matrix 
#' from a simulated single cell RNA sequencing data with two conditions
#'
#' @format A list with 5 items:
#' \describe{
#'   \item{numb_cond1}{the count matrix of condition 1}
#'   \item{numb_cond2}{the count matrix of condition 2}
#'   \item{knn_mat}{the similariy matrix}
#'   \item{knnGraphs}{the knn graphs information used to calculate the 
#'   similarity matrix.}
#'   \item{labels}{the clusters' label for each simulated single cell}
#' }
"simulation"

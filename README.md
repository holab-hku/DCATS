
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DCATS

<!-- badges: start -->
<!-- badges: end -->

This R package contains methods to detect the differential composition
abundances between multiple conditions in singel-cell experiments.

The **latest** version of the `DCATS` package is 0.99.0.

## Installation

### From Biocounductor

``` r
if (!requireNamespace("BiocManager"))
install.packages("BiocManager")
BiocManager::install("DCTAS")
```

### From R

The **latest** `DCATS` package can be conveniently installed using the
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/)
package thus:

``` r
## install dependencies
install.packages(c("MCMCpack", "matrixStats", "robustbase", "aod", "e1071"))
## dependencies for vignette
install.packages(c("SeuratObject", "Seurat", "robustbase", "aod", "e1071"))
devtools::install_github('satijalab/seurat-data')
```

``` r
# install.packages("devtools")
devtools::install_github("holab-hku/DCATS", build_vignettes = TRUE)
```

You can also install `DCATS` without building the vignette:

    devtools::install_github("holab-hku/DCATS")

#### For development

Download this repository to your local machine and open it in Rstudio as
a project, and build it by install and restart.

## Getting started

The best place to start are the vignettes. Inside an R session, load
`DCATS` and then browse the vignette about the usage guidance of
`DCATS`:

``` r
library(DCATS)
browseVignettes("DCATS")
```

The tutorial demonstrating how to use DCATS after using
[`Seurat`](https://satijalab.org/seurat/index.html) pipeline to process
data can be found in
[`Integrate DCATS with Seurat pipeline`](https://htmlpreview.github.io/?https://github.com/linxy29/DCATS_anlysis/blob/master/vignette/Integrate_with_seurat.html).

### Example

This is a basic example which shows you how to estimate a similarity
matrix from KNN graph and do the differential abundance test using this
similarity matrix.

``` r
library(DCATS)
data("simulation")
knn_mat = knn_simMat(simulation$knnGraphs, simulation$labels)
sim_count = rbind(simulation$numb_cond1, simulation$numb_cond2)
sim_design = data.frame(condition = c("c1", "c1", "c2"))
knn_mat[colnames(sim_count),]
res = dcats_GLM(as.matrix(sim_count), sim_design, similarity_mat = knn_mat)
print(res$LRT_pvals)
```

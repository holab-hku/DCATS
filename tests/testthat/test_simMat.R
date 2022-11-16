library(testthat)
library(DCATS)

test_that("create_simMat generate a matrix we want", {
  stdRes = matrix(c(0.9, rep(0.1/3, 4), 0.9, rep(0.1/3, 4), 0.9, rep(0.1/3, 4), 0.9), ncol = 4)
  
  expect_equal(create_simMat(4, 0.1), stdRes)
  
})


#> Test passed 🌈
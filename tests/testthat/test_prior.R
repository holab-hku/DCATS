library(testthat)
library(DCATS)

test_that("whether simulator works", {
  K <- 2
  totals1 = c(100, 800, 1300, 600)
  totals2 = c(250, 700, 1100)
  diri_s1 = diri_s2 = rep(1, K) * 20
  simil_mat = create_simMat(K, 0.2)
  sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  
  expect_equal(dim(sim_dat$numb_cond1)[2], K)
  expect_equal(dim(sim_dat$numb_cond2)[2], K)
  expect_equal(rowSums(sim_dat$numb_cond1), totals1)
  expect_equal(rowSums(sim_dat$numb_cond2), totals2)
})

test_that("whether EM process works", {
  X = c(100, 300, 1500, 500, 1000)
  simMM = create_simMat(5, confuse_rate=0.2)
  res = multinom_EM(X, simMM)
  
  expect_equal(sum(res$mu), 1)
})

test_that("whether phi estimation works", {
  set.seed(123)
  K <- 3
  totals1 = c(100, 800, 1300, 600)
  totals2 = c(250, 700, 1100)
  diri_s1 = rep(1, K) * 20
  diri_s2 = rep(1, K) * 20
  simil_mat = DCATS::create_simMat(K, confuse_rate=0.2)
  sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  sim_count = rbind(sim_dat$numb_cond1, sim_dat$numb_cond2)
  sim_design = data.frame(condition = c("g1", "g1", "g1", "g1", "g2", "g2", "g2"))
  phi = DCATS::getPhi(sim_count, sim_design)
  
  expect_equal(as.numeric(round(phi, 6)), 0.007930)
})

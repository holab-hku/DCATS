library(testthat)
library(DCATS)

test_that("whether dcats_GLM works as we expected", {
  set.seed(3)
  K <- 3
  totals1 = rep(300, 4)
  totals2 = rep(300, 3)
  diri_s1 = rep(1/K, K) * 200
  diri_s2 = c(1/3, 1/2, 1/6) * 200
  simil_mat = create_simMat(K, confuse_rate=0.1)
  sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  sim_count = rbind(sim_dat$numb_cond1, sim_dat$numb_cond2)
  sim_design = matrix(c("g1", "g1", "g1", "g1", "g2", "g2", "g2"), ncol = 1)
  res = dcats_GLM(sim_count, sim_design)

  expect_equal(as.vector(res$LRT_pvals) > 0.05, c(TRUE, FALSE, FALSE))

})


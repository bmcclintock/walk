# Tests for the joint (multi-individual) likelihood. Simulated data only.
# Thread use is capped at 2 (CRAN policy).

test_that("bvn_rect_prob (used by proc_Langevin) matches mvtnorm::pmvnorm", {
  skip_if_not_installed("mvtnorm")
  set.seed(3)
  l1 <- stats::runif(60, -3, 2); l2 <- stats::runif(60, -3, 2)
  for (r in c(-0.999, -0.5, 0, 0.7, 0.999)) {
    a <- walk:::bvn_rect_prob(l1, l1 + 1, l2, l2 + 1, rep(r, 60))
    S <- matrix(c(1, r, r, 1), 2)
    b <- vapply(seq_along(l1), function(i)
      mvtnorm::pmvnorm(c(l1[i], l2[i]), c(l1[i], l2[i]) + 1, sigma = S, keepAttr = FALSE), 0)
    expect_lt(max(abs(a - b)), 1e-7)
  }
})

test_that("proc_Langevin: dt in hours, within individual, with segments", {
  w <- sim_world()
  pA <- suppressWarnings(proc_Langevin(w$A, w$cv))
  expect_equal(attr(pA$times, "time_unit"), "hours")
  expect_equal(pA$times$dt, c(0, diff(as.numeric(w$A$date))) / 3600)
  pAB <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv))
  expect_true(all(pAB$times$dt >= 0))
  expect_equal(pAB$times$dt[nrow(w$A) + 1], 0)
  expect_equal(unique(pAB$times$seg), 1:2)
  expect_equal(unique(pAB$times$id), c("A", "B"))
  pg <- suppressWarnings(proc_Langevin(w$A, w$cv, max_gap = 2))
  expect_equal(max(pg$times$seg), sum(diff(as.numeric(w$A$date)) / 3600 > 2) + 1)
  expect_true(all(pg$times$dt <= 2))
  expect_error(proc_Langevin(as_dl(rbind(w$A[1:5, ], w$B[1:5, ], w$A[6:10, ])), w$cv),
               "contiguous")
})

test_that("joint likelihood: equals the original forward algorithm and is additive", {
  w <- sim_world(); par <- par_test()
  pA <- suppressWarnings(proc_Langevin(w$A, w$cv))
  pB <- suppressWarnings(proc_Langevin(w$B, w$cv))
  dA <- dl_of(make_walk_data(pA, w$cv))
  dB <- dl_of(make_walk_data(pB, w$cv))
  vA <- ctmc_n2ll_joint(par, dA)
  expect_equal(as.numeric(vA), ctmc_n2ll(par, dA), tolerance = 1e-10)
  vB <- ctmc_n2ll_joint(par, dB)
  pAB <- bind_proc_data(list(pA, pB))
  dAB <- dl_of(make_walk_data(pAB, w$cv))
  expect_equal(as.numeric(ctmc_n2ll_joint(par, dAB)), as.numeric(vA + vB), tolerance = 1e-10)
  # bind_proc_data() == one proc_Langevin() call on the stacked data
  p1 <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv))
  expect_equal(p1$L, pAB$L)
  expect_equal(p1$times$dt, pAB$times$dt)
})

test_that("1 and 2 threads give identical results", {
  w <- sim_world(); par <- par_test()
  pAB <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, max_gap = 3))
  d1 <- dl_of(make_walk_data(pAB, w$cv), nthreads = 1)
  d2 <- d1; d2$nthreads <- 2L
  expect_identical(ctmc_n2ll_joint(par, d1, by_segment = TRUE),
                   ctmc_n2ll_joint(par, d2, by_segment = TRUE))
  # asking for more threads than allowed is harmless
  d9 <- d1; d9$nthreads <- 64L
  expect_identical(as.numeric(ctmc_n2ll_joint(par, d9)), as.numeric(ctmc_n2ll_joint(par, d1)))
})

test_that("skip_zero treats impossible fixes as missing", {
  w <- sim_world(); par <- par_test()
  A <- w$A
  # Fix 15: a precise (GPS-like) fix 1 second after fix 14 but >= 15 km away, still
  # on the grid. With uniformisation tolerance 1e-8 the chain cannot cover that
  # distance in 1 s, so the fix has zero predictive probability.
  A$date[15] <- A$date[14] + 1
  A$x[15] <- if (A$x[14] < 30e3) A$x[14] + 18e3 else A$x[14] - 18e3
  A$y[15] <- A$y[14]
  A$smaj[15] <- NA; A$smin[15] <- NA; A$eor[15] <- NA; A$x.err[15] <- 50; A$y.err[15] <- 50
  expect_true(all(diff(as.numeric(A$date)) >= 0))
  wd <- make_walk_data(suppressWarnings(proc_Langevin(A, w$cv)), w$cv)
  d <- dl_of(wd)
  expect_equal(as.numeric(ctmc_n2ll_joint(par, d)), Inf)
  d$skip_zero <- TRUE
  v <- ctmc_n2ll_joint(par, d)
  expect_equal(attr(v, "zero_obs"), 15L)
  wd2 <- wd; wd2$L[15, ] <- 0; wd2$L <- Matrix::drop0(wd2$L)
  expect_equal(as.numeric(v), as.numeric(ctmc_n2ll_joint(par, dl_of(wd2))), tolerance = 1e-10)
})

test_that("joint fit converges and predict_ctmc returns distributions", {
  skip_on_cran()
  w <- sim_world()
  pAB <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv))
  wd <- add_d(make_walk_data(pAB, w$cv))
  f <- suppressMessages(suppressWarnings(fit_ctmc(wd, model_parameters = mp_test(), start = list(
    beta_q_r = par_test()[1], beta_q_m = c(0, 0)), nthreads = 2,
    control = list(eval.max = 200, iter.max = 100))))
  expect_equal(f$opt$convergence, 0)
  expect_equal(f$n_seg, 2L)
  expect_true(all(is.finite(f$vcov)))
  pr <- predict_ctmc(f, wd, segments = "B")
  expect_equal(unname(Matrix::rowSums(pr$local_state_prob)),
               rep(1, nrow(pr$local_state_prob)), tolerance = 1e-8)
})

test_that("reference ctmc_n2ll equals the joint likelihood when the first fix has no information", {
  w <- sim_world(); par <- par_test()
  wd <- make_walk_data(suppressWarnings(proc_Langevin(w$A, w$cv)), w$cv)
  wd$L[1, ] <- 0; wd$L <- Matrix::drop0(wd$L)
  d <- dl_of(wd)
  expect_true(is.finite(ctmc_n2ll(par, d)))
  expect_equal(ctmc_n2ll(par, d), as.numeric(ctmc_n2ll_joint(par, d)), tolerance = 1e-10)
})

# Audit 2: argument checks, consistency fixes and streamlining helpers.

w2 <- function() {
  w <- sim_world(); w$p0 <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv)); w
}
st2 <- function() list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1])

test_that("make_walk_data: clear error for data.frame proc output; dense L accepted", {
  w <- w2()
  pdf <- suppressWarnings(proc_Langevin(as_dl(w$A), w$cv, return_type = "data.frame"))
  expect_error(make_walk_data(pdf, w$cv), "proc_data")
  pd <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, return_type = "dense"))
  wd <- suppressWarnings(make_walk_data(pd, w$cv))
  expect_s4_class(wd$L, "dgCMatrix")
})

test_that("difference covariates are added automatically (add_diff)", {
  w <- w2()
  wd <- suppressWarnings(make_walk_data(w$p0, w$cv))
  expect_equal(wd$q_m$d_c1, wd$q_m$c1 - wd$q_m$from_c1)
  expect_equal(wd$q_m$d_c2, wd$q_m$c2 - wd$q_m$from_c2)
  expect_false("d_c1" %in% names(suppressWarnings(make_walk_data(w$p0, w$cv, add_diff = FALSE))$q_m))
  # usable directly, no manual step
  f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2())
  expect_true(is.finite(as.numeric(ctmc_n2ll_joint(f$par, f$data_list))))
})

test_that("gradients do not depend on clipRast windows", {
  w <- w2()
  pc <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, clipRast = TRUE))
  wf <- suppressWarnings(make_walk_data(w$p0, w$cv, grad = "c1"))
  for (wc in list(suppressWarnings(make_walk_data(pc, w$cv, grad = "c1")),
                  suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, clipRast = 6000)), w$cv, grad = "c1")))) {
    m <- match(paste(wc$q_m$from_cell, wc$q_m$cell), paste(wf$q_m$from_cell, wf$q_m$cell))
    expect_false(anyNA(m))
    expect_equal(wc$q_m$c1_grad, wf$q_m$c1_grad[m], tolerance = 1e-10)
  }
  expect_error(make_walk_data(w$p0, w$cv, grad = "nope"), "grad")
})

test_that("ctmc_model / ctmc_control validate arguments early", {
  expect_error(ctmc_model(link = "bogus"), "link")
  expect_error(ctmc_model(a = 0.5), "'a'")
  expect_error(ctmc_model("~1"), "formula")
  expect_error(ctmc_model(link = "logit", L = 2, U = 1), "U>L")
  expect_error(ctmc_control(form = "bogus"), "form")
  expect_error(ctmc_control(form = "add"), "not implemented")
  expect_error(ctmc_control(delta = "bogus"), "delta")
  expect_error(ctmc_control(q_r = list(form = ~1)), "ctmc_model")
  expect_silent(ctmc_control(q_r = ctmc_model(~1, link = "logit", L = 0, U = 5), delta = "stationary"))
})

test_that("fit_ctmc checks 'start'", {
  w <- w2(); wd <- suppressWarnings(make_walk_data(w$p0, w$cv))
  expect_error(fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = list(beta_q_r = 0, beta_q_m = 1)), "length 1")
  expect_error(fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = list(beta_r = 0)), "Unknown")
  expect_error(fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = c(0, 0, 0)), "named list")
  f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = list(beta_q_m = c(0.3, -0.2)))
  expect_equal(unname(f$par), c(0, 0.3, -0.2))
})

test_that("get_lim_ud hpd sets hold at least the nominal mass", {
  w <- w2(); wd <- suppressWarnings(make_walk_data(w$p0, w$cv))
  f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2())
  u <- get_lim_ud(f, hpd = c(0.5, 0.95))
  for (lv in c(50, 95)) {
    inset <- !is.na(u[[paste0("hpd_", lv)]])
    expect_gte(sum(u$ud[inset]), lv / 100 - 1e-12)
    # smallest such set: dropping its least likely cell goes below the level
    expect_lt(sum(u$ud[inset]) - min(u$ud[inset]), lv / 100)
  }
})

test_that("predict_ctmc: rows sum to 1 after trunc_tol; works without walk_data; aux warnings", {
  w <- w2(); wd <- suppressWarnings(make_walk_data(w$p0, w$cv))
  f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2())
  p1 <- predict_ctmc(f, wd, trunc_tol = 1e-3)
  expect_true(all(abs(Matrix::rowSums(p1$local_state_prob) - 1) < 1e-10))
  p2 <- predict_ctmc(f)            # no walk_data needed
  p3 <- predict_ctmc(f, wd)
  expect_equal(p2$local_state_prob, p3$local_state_prob)
  expect_equal(p2$times, p3$times)
  # aux before the first fix of every individual -> warning, not silent
  early <- min(wd$times$timestamp) - 86400
  expect_warning(predict_ctmc(f, aux_timestamp = early), "not predicted")
  # max_gap: aux time inside the gap
  pg <- suppressWarnings(proc_Langevin(as_dl(w$A), w$cv, max_gap = 2))
  wg <- suppressWarnings(make_walk_data(pg, w$cv))
  fg <- fit_ctmc(wg, model_parameters = mp_test(), fit = FALSE, start = st2())
  tt <- wg$times; k <- which(tt$seg[-1] != tt$seg[-nrow(tt)])[1]
  mid <- tt$timestamp[k] + (tt$timestamp[k + 1] - tt$timestamp[k]) / 2
  expect_warning(pr <- predict_ctmc(fg, aux_timestamp = mid), "max_gap")
  expect_equal(sum(pr$times$type == "p"), 0)
})

test_that("walk_rast maps UD and predictions to the user's raster (also when clipped)", {
  w <- w2()
  pc <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, clipRast = TRUE))
  for (wd in list(suppressWarnings(make_walk_data(w$p0, w$cv)), suppressWarnings(make_walk_data(pc, w$cv)))) {
    f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2())
    u <- get_lim_ud(f, hpd = 0.5)
    r <- walk_rast(u)
    expect_true(terra::compareGeom(r, w$cv))
    expect_equal(names(r), c("ud", "hpd_50"))
    expect_equal(as.vector(terra::values(r[[1]]))[u$cell], u$ud)
    expect_equal(sum(terra::values(r[[1]]), na.rm = TRUE), 1)
    pr <- predict_ctmc(f)
    rp <- walk_rast(pr, rows = 1:3)
    expect_equal(terra::nlyr(rp), 3)
    expect_equal(unname(terra::global(rp, "sum", na.rm = TRUE)[, 1]), rep(1, 3), tolerance = 1e-8)
  }
  expect_error(walk_rast(data.frame(a = 1)), "get_lim_ud")
})

test_that("unknown arguments are errors, not silently ignored", {
  w <- w2()
  # make_walk_data() and predict_ctmc() have no '...'
  expect_error(make_walk_data(w$p0, w$cv, bogus = 5000), "unused argument")
  wd <- suppressWarnings(make_walk_data(w$p0, w$cv))
  f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2())
  expect_error(predict_ctmc(f, foo = 1), "unused argument")
  # fit_ctmc() passes only optimx::optimr()'s own arguments through '...'
  expect_error(fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2(), contorl = list(maxit = 1)),
               "Unused argument\\(s\\) in fit_ctmc\\(\\): contorl")
  expect_silent(fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = st2(), control = list(maxit = 1)))
  # proc_*: 'return_type' must be valid
  expect_error(proc_Langevin(as_dl(w$A), w$cv, return_type = "hours"), "must be one of")
  expect_error(proc_Langevin(as_dl(w$A), w$cv, return_type = "sparce"), "must be one of")
  expect_error(proc_telem(mock_telemetry(w$A, w$cv), w$cv, return_type = 1), "must be one of")
  m <- suppressWarnings(proc_Langevin(as_dl(w$A), w$cv, return_type = "matrix_df"))
  expect_equal(colnames(m), c("obs", "cell", "lik"))
})

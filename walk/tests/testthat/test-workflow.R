# Whole workflow on simulated data:
# proc_Langevin/proc_telem -> make_walk_data -> ctmc_control -> fit_ctmc -> get_lim_ud -> predict_ctmc
# with clipping, several individuals, masks, and time units.

wf_world <- function() {
  w <- sim_world()
  xy <- terra::xyFromCell(w$cv, seq_len(terra::ncell(w$cv)))
  mask <- terra::rast(w$cv[[1]]); terra::values(mask) <- 1
  mask[xy[, 1] > 50e3 & xy[, 2] > 30e3] <- NA
  w$mask <- mask
  w$st <- list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1])
  w
}
wf_fit <- function(wd, ...) suppressMessages(fit_ctmc(add_d(wd), model_parameters = mp_test(), fit = FALSE,
                                                     start = list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1]), ...))

test_that("cell numbers refer to the supplied cell_data after clipRast (TRUE and buffer)", {
  w <- wf_world()
  AB <- as_dl(rbind(w$A, w$B))
  pf <- suppressWarnings(proc_Langevin(AB, w$cv))
  wf <- suppressWarnings(make_walk_data(pf, w$cv, rast_mask = w$mask))
  ff <- wf_fit(wf); uf <- get_lim_ud(ff)
  for (k in 1:3) {
    wd <- switch(k,
      suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(AB, w$cv, clipRast = TRUE)), w$cv, rast_mask = w$mask)),
      suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(AB, w$cv, clipRast = 6000)), w$cv, rast_mask = w$mask)),
      suppressWarnings(make_walk_data(suppressWarnings(proc_telem(list(mock_telemetry(w$A, w$cv), mock_telemetry(w$B, w$cv)),
                                                                  w$cv, clipRast = TRUE)), w$cv, rast_mask = w$mask)))
    expect_lt(nrow(wd$q_r), nrow(wf$q_r))
    cd <- terra::unwrap(wd$cell_data)
    expect_equal(terra::ncell(cd), terra::ncell(w$cv))
    # q_r / q_m cell numbers are cells of w$cv with the q_r coordinates and covariates
    expect_equal(terra::xyFromCell(w$cv, wd$q_r$cell), as.matrix(wd$q_r[, c("x", "y")]), ignore_attr = TRUE)
    expect_equal(terra::values(w$cv)[wd$q_r$cell, "c1"], wd$q_r$c1)
    expect_equal(terra::values(w$cv)[wd$q_m$from_cell, "c2"], wd$q_m$from_c2)
    expect_true(all(wd$q_m$cell %in% wd$q_r$cell))
    expect_true(all(wd$q_r$cell %in% wf$q_r$cell))
    # likelihood columns: L[, cellx] of the window equals the full-raster column of that cell
    if (k == 2) expect_equal(wd$L, wf$L[, match(wd$q_r$cell, wf$q_r$cell)])
    f <- wf_fit(wd); u <- get_lim_ud(f)
    expect_equal(u$cell, wd$q_r$cell[u$cellx])
    # mapping onto the raster the user has puts the UD on the right cells
    r <- terra::rast(w$cv[[1]]); r[] <- NA; r[u$cell] <- u$ud
    pk <- u$cellx[which.max(u$ud)]
    expect_equal(terra::xyFromCell(r, which.max(terra::values(r))), as.matrix(wd$q_r[pk, c("x", "y")]), ignore_attr = TRUE)
    pr <- predict_ctmc(f, wd)
    expect_identical(pr$cell_map, f$data_list$cell_map)
    expect_equal(ncol(pr$local_state_prob), nrow(pr$cell_map))
  }
})

test_that("NA cells and NA covariates in cell_data", {
  w <- wf_world()
  p <- suppressWarnings(proc_Langevin(w$A, w$cv))
  i <- terra::cellFromXY(w$cv, cbind(45e3, 10e3))
  # NA in the first layer: inadmissible (as in proc_Langevin), dropped from q_r
  cv1 <- w$cv; cv1[[1]][i] <- NA
  wd <- suppressWarnings(make_walk_data(p, cv1))
  expect_false(i %in% wd$q_r$cell); expect_equal(nrow(wd$q_r), terra::ncell(w$cv) - 1)
  expect_false(anyNA(wd$q_r$c1))
  f <- wf_fit(wd); expect_equal(nrow(f$data_list$X_q_m), nrow(wd$q_m))
  # NA in another layer used by the model: warning in make_walk_data, error (not silent row loss) in fit_ctmc
  cv2 <- w$cv; cv2[[2]][i] <- NA
  expect_warning(wd2 <- make_walk_data(p, cv2), "NA values")
  expect_error(wf_fit(wd2), "NA")
  # masking the cell fixes it
  m <- w$mask; m[i] <- NA
  wd3 <- suppressWarnings(make_walk_data(p, cv2, rast_mask = m))
  expect_equal(nrow(wf_fit(wd3)$data_list$X_q_m), nrow(wd3$q_m))
})

test_that("delta options in fit_ctmc", {
  w <- wf_world()
  wd <- suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv)), w$cv, rast_mask = w$mask))
  par <- par_test()
  fu <- wf_fit(wd)
  mpn <- mp_test(); mpn$delta <- rep(2, nrow(wd$q_r))
  fn <- suppressMessages(fit_ctmc(add_d(wd), model_parameters = mpn, fit = FALSE))
  expect_equal(fn$data_list$delta, fu$data_list$delta)
  expect_equal(ctmc_n2ll_joint(par, fn$data_list), ctmc_n2ll_joint(par, fu$data_list))
  mpn$delta <- rep(1, 3); expect_error(fit_ctmc(add_d(wd), model_parameters = mpn, fit = FALSE), "number of cells")
  mpn$delta <- "unifrom"; expect_error(fit_ctmc(add_d(wd), model_parameters = mpn, fit = FALSE), "delta")
  mps <- mp_test(); mps$delta <- "stationary"
  fs <- suppressMessages(fit_ctmc(add_d(wd), model_parameters = mps, fit = FALSE))
  expect_true(is.finite(ctmc_n2ll_joint(par, fs$data_list)))
})

test_that("fit_ctmc segment summary and get_lim_ud clip", {
  w <- wf_world()
  wd <- suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, max_gap = 3)), w$cv))
  f <- suppressMessages(fit_ctmc(add_d(wd), model_parameters = mp_test(), hessian = FALSE,
                                 start = list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1]),
                                 control = list(iter.max = 2, eval.max = 3)))
  expect_equal(nrow(f$seg_info), f$n_seg)
  expect_equal(f$seg_info$id, wd$times$id[f$seg_info$start])
  expect_equal(sum(f$seg_info$n), nrow(wd$times))
  expect_equal(names(f$seg_n2ll), as.character(f$seg_info$seg))
  expect_equal(sum(f$seg_n2ll), ctmc_n2ll_joint(f$par, f$data_list))
  # clip argument of get_lim_ud is used
  mpc <- mp_test(); mpc$clip <- 0.05
  fc <- suppressMessages(fit_ctmc(add_d(wd), model_parameters = mpc, fit = FALSE,
                                  start = list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1])))
  f0 <- wf_fit(wd)
  expect_false(isTRUE(all.equal(get_lim_ud(fc, clip = TRUE)$ud, get_lim_ud(fc, clip = FALSE)$ud)))
  expect_equal(get_lim_ud(fc, clip = FALSE)$ud, get_lim_ud(f0)$ud)
})

test_that("predict_ctmc: segments, aux_timestamp, and walk_data checks", {
  w <- wf_world()
  p <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv))
  wd <- suppressWarnings(make_walk_data(p, w$cv, rast_mask = w$mask))
  f <- wf_fit(wd)
  # factor ids
  expect_equal(unique(predict_ctmc(f, wd, segments = factor("B"))$times$id), "B")
  expect_warning(predict_ctmc(f, wd, segments = c("B", "Z")), "not found")
  expect_error(suppressWarnings(predict_ctmc(f, wd, segments = 9)), "No segments")
  # aux time equal to an observation time is not duplicated, and obs dt is unchanged
  tA <- wd$times[wd$times$id == "A", ]
  pr <- predict_ctmc(f, wd, aux_timestamp = c(tA$timestamp[5], tA$timestamp[5] + 600), segments = "A")
  expect_equal(nrow(pr$times), nrow(tA) + 1)
  expect_false(any(duplicated(pr$times$timestamp)))
  expect_equal(pr$times$dt[pr$times$type == "o"][-6], tA$dt[-6], tolerance = 1e-10)
  # numeric (seconds since epoch) and other time zones are converted to the data time stamps
  a <- tA$timestamp[3] + 600
  r1 <- predict_ctmc(f, wd, aux_timestamp = a, segments = "A")
  r2 <- predict_ctmc(f, wd, aux_timestamp = as.numeric(a), segments = "A")
  a3 <- a; attr(a3, "tzone") <- "America/Los_Angeles"
  r3 <- expect_silent(predict_ctmc(f, wd, aux_timestamp = a3, segments = "A"))
  r4 <- predict_ctmc(f, wd, aux_timestamp = data.frame(id = factor("A"), timestamp = as.POSIXlt(a)))
  for (r in list(r2, r3)) {
    expect_equal(r$times, r1$times); expect_equal(r$local_state_prob, r1$local_state_prob)
  }
  expect_equal(attr(r3$times$timestamp, "tzone"), attr(tA$timestamp, "tzone"))
  expect_equal(sum(r4$times$type == "p"), 1)
  # walk_data that does not belong to the fit
  wc <- suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, clipRast = TRUE)),
                                        w$cv, rast_mask = w$mask))
  expect_error(predict_ctmc(f, wc), "does not match")
  wm <- suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, time_unit = "mins")),
                                        w$cv, rast_mask = w$mask))
  expect_error(predict_ctmc(f, wm), "time_unit")
})

test_that("make_subset_neigh treats each segment start like the first observation", {
  w <- wf_world()
  wd <- suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv)), w$cv))
  s <- make_subset_neigh(wd, vmax = 2000)
  N <- nrow(wd$times); first <- which(c(TRUE, diff(wd$times$seg) != 0))
  len <- diff(s$active_offsets)
  expect_equal(length(len), N)
  expect_equal(len[first], rep(1, length(first)))
  expect_true(all(len[-first] > 1))
  # a single segment gives the same sets as that animal alone
  sA <- make_subset_neigh(suppressWarnings(make_walk_data(suppressWarnings(proc_Langevin(w$A, w$cv)), w$cv)), vmax = 2000)
  nA <- sum(wd$times$id == "A")
  expect_equal(s$active_indices[seq_len(s$active_offsets[nA + 1])], sA$active_indices)
})

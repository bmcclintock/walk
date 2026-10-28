# proc_telem() / proc_Langevin() / bind_proc_data() consistency.
# ctmm telemetry objects are mocked (helper-sim.R: mock_telemetry) so ctmm is not needed.

test_that("proc_telem and proc_Langevin have the same arguments", {
  a <- names(formals(proc_telem)); b <- names(formals(proc_Langevin))
  expect_identical(a, b)
  # same defaults except time_unit: "hours" for proc_telem, NULL (the data's unit) for proc_Langevin
  ft <- formals(proc_telem)[-1]; fl <- formals(proc_Langevin)[-1]
  expect_identical(ft[names(ft) != "time_unit"], fl[names(fl) != "time_unit"])
  expect_identical(ft$time_unit, "hours"); expect_null(fl$time_unit)
})

test_that("proc_telem: id from ctmm identity, list input, max_gap, same times as proc_Langevin", {
  w <- sim_world()
  tA <- mock_telemetry(w$A, w$cv); tB <- mock_telemetry(w$B, w$cv)
  pA <- suppressWarnings(proc_telem(tA, w$cv))
  expect_named(pA$times, c("obs", "id", "seg", "timestamp", "dt"))
  expect_equal(unique(pA$times$id), "A")                   # from attr(, "info")$identity
  expect_equal(attr(pA$times, "time_unit"), "hours")
  expect_equal(pA$times$dt, c(0, diff(as.numeric(w$A$date))) / 3600)
  # list of telemetry objects (what ctmm::as.telemetry() returns for several animals)
  pAB <- suppressWarnings(proc_telem(list(tA, tB), w$cv, max_gap = 3))
  lAB <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, max_gap = 3))
  expect_identical(pAB$times[, c("obs", "id", "seg", "dt")], lAB$times[, c("obs", "id", "seg", "dt")])
  expect_equal(as.numeric(pAB$times$timestamp), as.numeric(lAB$times$timestamp))
  expect_true(all(pAB$times$dt <= 3))
  # same locations and error covariance -> same cell likelihoods (up to lon/lat round trip)
  expect_equal(pAB$L, lAB$L, tolerance = 1e-6)
  # identity missing -> list names, then position; single object -> "1"; 'id' column wins
  tN <- tA; attr(tN, "info")$identity <- NULL
  expect_equal(unique(suppressWarnings(proc_telem(list(foo = tN), w$cv))$times$id), "foo")
  expect_equal(unique(suppressWarnings(proc_telem(list(tN), w$cv))$times$id), "1")
  expect_equal(unique(suppressWarnings(proc_telem(tN, w$cv))$times$id), "1")
  tC <- tA; tC$id <- "zz"
  expect_equal(unique(suppressWarnings(proc_telem(tC, w$cv))$times$id), "zz")
  expect_equal(unique(suppressWarnings(proc_telem(tC, w$cv, id = NULL))$times$id), "A")
  # the same animal twice is not a valid list (non-contiguous after stacking is impossible,
  # but an id repeated in separate elements is)
  expect_error(proc_telem(list(tA, tB, tA), w$cv), "contiguous")
  # an invalid return_type is rejected
  expect_error(proc_telem(tA, w$cv, "days"), "must be one of")
  expect_error(proc_telem(w$A, w$cv), "telemetry")
})

test_that("proc_Langevin accepts a list of tracks", {
  w <- sim_world()
  l1 <- suppressWarnings(proc_Langevin(list(w$A, w$B), w$cv, max_gap = 3))
  l2 <- suppressWarnings(proc_Langevin(as_dl(rbind(w$A, w$B)), w$cv, max_gap = 3))
  expect_identical(l1$times, l2$times)
  expect_equal(l1$L, l2$L)
  A <- w$A; A$id <- NULL
  expect_equal(unique(suppressWarnings(proc_Langevin(list(a1 = as_dl(A)), w$cv))$times$id), "a1")
})

test_that("bind_proc_data combines proc_telem, proc_Langevin, or a mix", {
  w <- sim_world(); par <- par_test()
  C <- sim_track("C", 20, 13, c(30e3, 18e3))
  tA <- mock_telemetry(w$A, w$cv); tB <- mock_telemetry(w$B, w$cv)
  ex <- terra::ext(w$cv)
  ex2 <- terra::ext(ex[1] + 2 * terra::res(w$cv)[1], ex[2], ex[3], ex[4] - 2 * terra::res(w$cv)[2])
  for (clip in list(FALSE, ex2)) {
    pA <- suppressWarnings(proc_telem(tA, w$cv, clipRast = clip, max_gap = 3))
    pB <- suppressWarnings(proc_telem(tB, w$cv, clipRast = clip, max_gap = 3))
    pC <- suppressWarnings(proc_Langevin(C, w$cv, clipRast = clip, max_gap = 3))
    # telem only == one proc_telem() call on the list
    b1 <- bind_proc_data(pA, pB)
    one <- suppressWarnings(proc_telem(list(tA, tB), w$cv, clipRast = clip, max_gap = 3))
    expect_identical(b1$times, one$times)
    expect_equal(b1$L, one$L)
    expect_identical(b1$clip, one$clip)
    # mix: telemetry + Langevin
    m <- bind_proc_data(list(pA, pB, pC))
    expect_equal(unique(m$times$id), c("A", "B", "C"))
    expect_equal(m$times$seg, c(pA$times$seg, max(pA$times$seg) + pB$times$seg,
                                max(pA$times$seg) + max(pB$times$seg) + pC$times$seg))
    expect_equal(m$times$obs, seq_len(nrow(m$L)))
    expect_identical(m$clip, pA$clip)
    cv <- if (isFALSE(clip)) w$cv else terra::crop(w$cv, terra::ext(m$clip$ext))
    # joint -2ll of the mix = sum of the parts
    f <- function(p) as.numeric(ctmc_n2ll_joint(par, dl_of(make_walk_data(p, cv))))
    expect_equal(f(m), f(pA) + f(pB) + f(pC), tolerance = 1e-10)
  }
  # mismatches are rejected
  pA <- suppressWarnings(proc_telem(tA, w$cv))
  pCd <- suppressWarnings(proc_Langevin(C, w$cv, time_unit = "days"))
  expect_error(bind_proc_data(pA, pCd), "time_unit")
  expect_error(bind_proc_data(pA, pA), "more than one input")
  pAc <- suppressWarnings(proc_telem(tA, w$cv, clipRast = TRUE))
  pCc <- suppressWarnings(proc_Langevin(C, w$cv, clipRast = TRUE))
  # different clipRast windows are combined on their bounding box (tested in test-domain.R)
  expect_false(identical(pAc$clip$ext, pCc$clip$ext))
  expect_silent(bind_proc_data(pAc, pCc))
  expect_error(bind_proc_data(pA, suppressWarnings(proc_Langevin(C, w$cv, return_type = "data.frame"))),
               "return_type")
  # 'ids' overrides
  b <- bind_proc_data(pA, suppressWarnings(proc_Langevin(C, w$cv)), ids = c("x", "y"))
  expect_equal(unique(b$times$id), c("x", "y"))
  expect_equal(unique(b$times$seg), 1:2)
})

test_that("make_walk_data accepts the full cell_data and rast_mask for clipped proc_data", {
  w <- sim_world()
  mask <- terra::rast(w$cv[[1]]); terra::values(mask) <- 1
  mask[terra::xyFromCell(mask, 1:terra::ncell(mask))[, 1] > 50e3] <- NA
  for (clip in list(TRUE, terra::ext(4e3, 46e3, 4e3, 36e3))) {
    p <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = clip))
    expect_lt(ncol(p$L), terra::ncell(w$cv))
    cvc <- terra::crop(w$cv, terra::ext(p$clip$ext)); mc <- terra::crop(mask, terra::ext(p$clip$ext))
    ref <- suppressWarnings(make_walk_data(p, cvc, rast_mask = mc))
    for (args in list(list(w$cv, mask), list(w$cv, mc), list(cvc, mask), list(w$cv, NULL))) {
      wd <- suppressWarnings(make_walk_data(p, args[[1]], rast_mask = args[[2]]))
      ref_k <- if (is.null(args[[2]])) suppressWarnings(make_walk_data(p, cvc)) else ref
      expect_identical(wd$L, ref_k$L)
      # same cells and model; 'cell' numbers refer to the raster supplied to make_walk_data()
      nc <- c("cell", "from_cell")
      expect_identical(wd$q_r[, setdiff(names(wd$q_r), nc)], ref_k$q_r[, setdiff(names(ref_k$q_r), nc)])
      expect_identical(wd$q_m[, setdiff(names(wd$q_m), nc)], ref_k$q_m[, setdiff(names(ref_k$q_m), nc)])
      cd <- terra::unwrap(wd$cell_data)
      expect_equal(terra::ncell(cd), terra::ncell(args[[1]]))
      expect_equal(terra::xyFromCell(cd, wd$q_r$cell), as.matrix(wd$q_r[, c("x", "y")]), ignore_attr = TRUE)
      expect_equal(terra::xyFromCell(cvc, ref_k$q_r$cell), as.matrix(ref_k$q_r[, c("x", "y")]), ignore_attr = TRUE)
    }
  }
  expect_error(make_walk_data(p, terra::aggregate(w$cv, 2)), "clipRast")
})

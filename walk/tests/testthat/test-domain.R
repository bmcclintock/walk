# Model domain: clipRast (TRUE / buffer / extent), bind_proc_data with different windows.

test_that("clipRast = <buffer> expands the error-neighbourhood box; TRUE equals 0", {
  w <- sim_world()
  p0 <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = TRUE))
  pz <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = 0))
  expect_identical(p0$clip, pz$clip)
  expect_equal(p0$L, pz$L)
  expect_equal(unname(p0$clip$full_dim), dim(w$cv)[1:2])
  p6 <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = 6000))
  e0 <- p0$clip$ext; e6 <- p6$clip$ext; re <- as.vector(terra::ext(w$cv)); rs <- terra::res(w$cv)[1]
  # every side moves out by the buffer (to the next cell edge), unless it hits the raster edge
  expect_true(all(e6[c(1, 3)] <= pmax(e0[c(1, 3)] - 6000 + 1e-6, re[c(1, 3)])))
  expect_true(all(e6[c(1, 3)] >  pmax(e0[c(1, 3)] - 6000 - rs, re[c(1, 3)]) - 1e-6))
  expect_true(all(e6[c(2, 4)] >= pmin(e0[c(2, 4)] + 6000 - 1e-6, re[c(2, 4)])))
  expect_gt(ncol(p6$L), ncol(p0$L))
  # the likelihood values are the same (no error neighbourhood is cut either way)
  expect_equal(sum(p6$L), sum(p0$L))
  expect_equal(sort(p6$L@x), sort(p0$L@x))
  # same as the equivalent extent
  pe <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = terra::ext(p6$clip$ext)))
  expect_identical(pe$clip, p6$clip); expect_equal(pe$L, p6$L)
  # a huge buffer = the whole raster (still a 'clip' element, but covering every cell)
  pb <- suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = 1e9))
  expect_equal(ncol(pb$L), terra::ncell(w$cv))
  # proc_telem too
  t6 <- suppressWarnings(proc_telem(mock_telemetry(w$A, w$cv), w$cv, clipRast = 6000))
  expect_equal(length(t6$clip$cell), ncol(t6$L))
  # invalid values
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "a", NA))
    expect_error(suppressWarnings(proc_Langevin(w$A, w$cv, clipRast = bad)), "clipRast")
})

test_that("bind_proc_data combines different clipRast windows on their bounding box", {
  w <- sim_world(); par <- par_test()
  A <- w$A; B <- w$B
  pA <- suppressWarnings(proc_Langevin(A, w$cv, clipRast = 4000))
  pB <- suppressWarnings(proc_Langevin(B, w$cv, clipRast = 4000))
  expect_false(identical(pA$clip, pB$clip))
  b <- bind_proc_data(pA, pB)
  ea <- pA$clip$ext; eb <- pB$clip$ext
  expect_equal(unname(b$clip$ext), unname(c(min(ea[1], eb[1]), max(ea[2], eb[2]), min(ea[3], eb[3]), max(ea[4], eb[4]))))
  expect_equal(ncol(b$L), length(b$clip$cell))
  # identical to processing each individual on the union window
  ex <- terra::ext(b$clip$ext)
  qA <- suppressWarnings(proc_Langevin(A, w$cv, clipRast = ex))
  qB <- suppressWarnings(proc_Langevin(B, w$cv, clipRast = ex))
  ref <- bind_proc_data(qA, qB)
  expect_identical(b$clip, ref$clip)
  expect_equal(b$L, ref$L)
  expect_identical(b$times, ref$times)
  # and the walk data / likelihood are the same
  d1 <- dl_of(make_walk_data(b, w$cv)); d2 <- dl_of(make_walk_data(ref, w$cv))
  expect_equal(as.numeric(ctmc_n2ll_joint(par, d1)), as.numeric(ctmc_n2ll_joint(par, d2)))
  # mixing a clipped and an unclipped input gives an unclipped result on the whole raster
  pBf <- suppressWarnings(proc_Langevin(B, w$cv))
  m <- bind_proc_data(pA, pBf)
  expect_null(m$clip)
  expect_equal(ncol(m$L), terra::ncell(w$cv))
  pAf <- suppressWarnings(proc_Langevin(A, w$cv))
  expect_equal(m$L, bind_proc_data(pAf, pBf)$L)
  # different rasters are still rejected
  cv2 <- terra::aggregate(w$cv, 2)
  expect_error(bind_proc_data(pA, suppressWarnings(proc_Langevin(B, cv2, clipRast = 4000))), "cell_data")
  expect_error(bind_proc_data(pA, suppressWarnings(proc_Langevin(B, cv2))), "cell_data")
  # a malformed clip element is rejected
  pA0 <- pA; pA0$clip$full_dim <- NULL
  expect_error(bind_proc_data(pA0, pB), "full_dim")
})

test_that("the domain only shifts the uniform-delta log-likelihood by log(n cells) per segment (mult, norm = FALSE)", {
  w <- sim_world(); par <- par_test()
  AB <- as_dl(rbind(w$A, w$B))
  pf <- suppressWarnings(proc_Langevin(AB, w$cv))
  pc <- suppressWarnings(proc_Langevin(AB, w$cv, clipRast = 4000))
  df <- dl_of(make_walk_data(pf, w$cv)); dc <- dl_of(make_walk_data(pc, w$cv))
  nseg <- length(unique(pf$times$seg))
  expect_lt(dc$ns, df$ns)
  # -2ll differs by 2 * log(n cells) per segment (the uniform delta); the walk does not reach
  # the edge of the 4 km buffer at these parameters, so nothing else changes
  d <- as.numeric(ctmc_n2ll_joint(par, df)) - as.numeric(ctmc_n2ll_joint(par, dc))
  expect_equal(d, 2 * nseg * log(df$ns / dc$ns), tolerance = 1e-6)
})

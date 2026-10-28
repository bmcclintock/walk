# Excluded observations: make_walk_data(exclude =, vmax =), exclude_obs(), and how fit_ctmc() /
# predict_ctmc() treat them.

outlier_world <- function(k = 12) {
  w <- sim_world(); A <- w$A
  A$x[k] <- A$x[k - 1] + 30e3 * sign(30e3 - A$x[k - 1]); A$y[k] <- A$y[k - 1]   # 30 km jump in 6 min
  A$x.err[k] <- A$y.err[k] <- 50; A$smaj[k] <- A$smin[k] <- A$eor[k] <- NA
  A$date[k] <- A$date[k - 1] + 360
  w$A <- A; w$k <- k
  w$p <- suppressWarnings(proc_Langevin(as_dl(rbind(A, w$B)), w$cv))
  w$wd <- add_d(suppressWarnings(make_walk_data(w$p, w$cv)))
  w
}
st0 <- function() list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1])
n2 <- function(wd, b = par_test(), delta = "uniform") {
  mp <- mp_test(); mp$delta <- delta
  as.numeric(ctmc_n2ll_joint(b, fit_ctmc(add_d(wd), model_parameters = mp, fit = FALSE)$data_list))
}
# delete row k and measure the next step from the previous retained row
delete_row <- function(wd, k) {
  tt <- wd$times; keep <- setdiff(seq_len(nrow(tt)), k)
  if (k < nrow(tt) && tt$seg[k + 1] == tt$seg[k]) tt$dt[k + 1] <- tt$dt[k] + tt$dt[k + 1]
  wd$L <- wd$L[keep, ]; wd$times <- tt[keep, ]; wd
}
pars <- list(par_test(), c(0, 0.5, 0.1), c(-3, 0, 0))

test_that("make_walk_data always returns times$excluded; exclude= and exclude_obs() agree", {
  w <- outlier_world()
  expect_true(all(is.na(w$wd$times$excluded)))
  expect_type(w$wd$times$excluded, "character")
  we <- add_d(suppressWarnings(make_walk_data(w$p, w$cv, exclude = w$k)))
  expect_equal(we$times$excluded[w$k], "user"); expect_true(all(is.na(we$times$excluded[-w$k])))
  expect_equal(sum(we$L[w$k, ]), 0); expect_equal(we$L[-w$k, ], w$wd$L[-w$k, ])
  wx <- exclude_obs(w$wd, w$k)
  expect_identical(wx$L, we$L); expect_identical(wx$times, we$times)
  # logical exclude, idempotent, reason, bad input
  expect_identical(suppressWarnings(make_walk_data(w$p, w$cv, exclude = seq_len(nrow(w$p$L)) == w$k))$L, we$L)
  expect_identical(exclude_obs(we, w$k, reason = "other")$times$excluded[w$k], "user")
  expect_equal(exclude_obs(w$wd, w$k, reason = "outlier")$times$excluded[w$k], "outlier")
  expect_error(make_walk_data(w$p, w$cv, exclude = 0), "row numbers")
  expect_error(make_walk_data(w$p, w$cv, exclude = c(TRUE, FALSE)), "one non-missing value")
  expect_error(make_walk_data(w$p, w$cv, vmax = -1), "vmax")
  expect_error(exclude_obs(w$wd, nrow(w$wd$L) + 1), "row numbers")
  expect_error(exclude_obs(w$wd, vmax = c(1, 2)), "vmax")
})

test_that("an excluded observation = deleting it, for interior and first observations", {
  w <- outlier_world(); wd <- w$wd
  first <- which(!duplicated(wd$times$seg))
  for (k in c(w$k, first)) {
    we <- exclude_obs(wd, k); wdl <- delete_row(wd, k)
    for (delta in c("uniform", "stationary")) for (b in pars)
      expect_equal(n2(we, b, delta), n2(wdl, b, delta), tolerance = 1e-9)
  }
  # several excluded rows at the start of a segment
  k <- first[2] + 0:2
  wdl <- wd; for (j in rev(k)) wdl <- delete_row(wdl, j)
  expect_equal(n2(exclude_obs(wd, k)), n2(wdl), tolerance = 1e-9)
  fe <- fit_ctmc(exclude_obs(wd, c(w$k, k)), model_parameters = mp_test(), fit = FALSE, start = st0())
  expect_equal(fe$seg_info$n_excluded, c(1, 3))
})

test_that("observations without cells are excluded as no_cells, exactly like exclude=", {
  w <- outlier_world()
  pA <- suppressWarnings(proc_Langevin(w$A, w$cv))
  m <- terra::rast(w$cv[[1]]); terra::values(m) <- 1; m[which(pA$L[5, ] > 0)] <- NA
  expect_warning(wm <- make_walk_data(pA, w$cv, rast_mask = m), "excluded.*A:5")
  nc <- which(wm$times$excluded == "no_cells")
  expect_true(5 %in% nc)
  expect_equal(unname(which(Matrix::rowSums(wm$L) == 0)), nc)
  # the same rows excluded by hand on the same (masked) domain give the same data and likelihood
  pm <- pA; pm$L[nc, ] <- 0; pm$L <- Matrix::drop0(pm$L)
  wu <- suppressWarnings(make_walk_data(pm, w$cv, rast_mask = m))
  expect_identical(wu$L, wm$L)
  wu2 <- make_walk_data(pA, w$cv, rast_mask = m, exclude = nc) |> suppressWarnings()
  expect_identical(wu2$L, wm$L); expect_equal(wu2$times$excluded[nc], rep("no_cells", length(nc)))
  expect_equal(n2(wm), n2(wu2))
})

test_that("vmax finds the speed outlier; make_walk_data(vmax=) == exclude_obs(vmax=)", {
  w <- outlier_world()
  s <- suppressWarnings(make_walk_data(w$p, w$cv, vmax = 15000))      # 15 km/h; the outlier needs ~300 km/h
  expect_equal(which(!is.na(s$times$excluded)), w$k)
  expect_equal(s$times$excluded[w$k], "speed")
  s2 <- exclude_obs(w$wd, vmax = 15000)
  expect_identical(s2$L, s$L); expect_identical(s2$times$excluded, s$times$excluded)
  expect_true(is.finite(n2(s)))
  # exclude= first, then vmax on the retained observations
  s3 <- suppressWarnings(make_walk_data(w$p, w$cv, exclude = w$k, vmax = 15000))
  expect_equal(s3$times$excluded[w$k], "user"); expect_equal(sum(!is.na(s3$times$excluded)), 1)
  expect_identical(exclude_obs(w$wd)$L, w$wd$L)  # nothing to do
})

test_that("an impossible observation stops fit_ctmc at the start values with a useful message", {
  w <- outlier_world()
  v <- ctmc_n2ll_joint(par_test(), fit_ctmc(w$wd, model_parameters = mp_test(), fit = FALSE)$data_list)
  expect_equal(as.numeric(v), Inf); expect_equal(attr(v, "zero_obs"), w$k)
  expect_error(suppressMessages(fit_ctmc(w$wd, model_parameters = mp_test(), start = st0())),
               paste0("likelihood is zero at the starting values: observation\\(s\\) ", w$k, " \\(A\\).*vmax"))
  # predict_ctmc on a fit = FALSE object at impossible parameters: same message
  f0 <- fit_ctmc(w$wd, model_parameters = mp_test(), fit = FALSE, start = st0())
  expect_error(predict_ctmc(f0), "likelihood is zero at the fitted parameters.*exclude_obs")
  # after excluding it everything works
  fe <- fit_ctmc(exclude_obs(w$wd, w$k), model_parameters = mp_test(), fit = FALSE, start = st0())
  pr <- predict_ctmc(fe)
  expect_false(anyNA(pr$local_state_prob@x))
  expect_equal(pr$times$excluded[pr$times$row_obs %in% w$k], "user")
  expect_equal(sum(!is.na(pr$times$excluded)), 1)
  rs <- Matrix::rowSums(pr$local_state_prob); expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("predict_ctmc: excluded first observations get the first retained observation's distribution", {
  w <- outlier_world(); first <- which(!duplicated(w$wd$times$seg))
  we <- exclude_obs(w$wd, c(w$k, first[2] + 0:1))
  fe <- fit_ctmc(we, model_parameters = mp_test(), fit = FALSE, start = st0())
  pr <- predict_ctmc(fe)
  r <- match(first[2] + 0:2, pr$times$row_obs)
  P <- pr$local_state_prob
  expect_equal(P[r[1], ], P[r[3], ], tolerance = 1e-10); expect_equal(P[r[2], ], P[r[3], ], tolerance = 1e-10)
  expect_equal(pr$times$excluded[r], c("user", "user", NA))
})

# ---- compiled vmax filter (audit 7) -----------------------------------------------------------
# random walk_ddl-like objects: several segments (non-contiguous labels), discs on a grid, outliers,
# zero time steps, empty rows, pre-excluded rows, and exact ties (integer coordinates / times)
rand_wd <- function(seed, N = 300, nx = 40, ny = 30) {
  set.seed(seed)
  seg <- sample(c(7, 2, 9, 4), N, replace = TRUE); seg <- seg[order(match(seg, unique(seg)))]
  seg <- rev(seg)
  dt <- sample(c(0, 1, 2, 3), N, replace = TRUE, prob = c(.1, .4, .3, .2))
  ncell <- nx * ny; cx <- sample(nx, N, replace = TRUE); cy <- sample(ny, N, replace = TRUE)
  ii <- jj <- vv <- list()
  for (k in seq_len(N)) {
    r <- sample(0:2, 1); o <- expand.grid(dx = -r:r, dy = -r:r)
    a <- cx[k] + o$dx; b <- cy[k] + o$dy; ok <- a >= 1 & a <= nx & b >= 1 & b <= ny
    ii[[k]] <- rep(k, sum(ok)); jj[[k]] <- (b[ok] - 1) * nx + a[ok]; vv[[k]] <- rep(1 / sum(ok), sum(ok))
  }
  L <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(vv), dims = c(N, ncell))
  L[sample(N, 5), ] <- 0                                 # empty rows
  L <- methods::as(L, "dgCMatrix")
  ex <- rep(NA_character_, N); ex[sample(N, 10)] <- "user"
  wd <- list(L = L, times = data.frame(seg = seg, dt = dt, excluded = ex),
             q_r = data.frame(x = rep(seq_len(nx), ny), y = rep(seq_len(ny), each = nx)),
             q_m = data.frame(h = rep(1, 4)))
  wd
}

test_that("compiled speed filter == plain-R reference (same rows, same order)", {
  ns <- asNamespace("walk")
  for (s in 1:25) {
    wd <- rand_wd(s)
    for (v in c(0.5, 1, 2, 5)) {
      a <- as.integer(ns$speed_filter_r(wd, v)); b <- ns$speed_filter(wd, v)
      expect_identical(b, a, info = paste("seed", s, "vmax", v))
    }
  }
  expect_gt(length(ns$speed_filter(rand_wd(1), 1)), 0)        # the comparison is not vacuous
  # no segment column: one segment; no q_m$h: zero buffer
  wd <- rand_wd(3); wd$times$seg <- NULL; wd$q_m <- data.frame(h = numeric(0))
  expect_identical(ns$speed_filter(wd, 1), as.integer(ns$speed_filter_r(wd, 1)))
})

test_that("row_support_bbox and exclude_rows match the Matrix-based versions", {
  ns <- asNamespace("walk")
  wd <- rand_wd(11); L <- wd$L; L@x[1:3] <- 0                  # stored zeros are ignored
  bb <- ns$row_support_bbox(L@p, L@i, L@x, nrow(L), wd$q_r$x, wd$q_r$y)
  T <- Matrix::summary(methods::as(L, "TsparseMatrix")); T <- T[T$x > 0, ]
  ref <- cbind(ns$tapply_full(wd$q_r$x[T$j], T$i, min, nrow(L)), ns$tapply_full(wd$q_r$x[T$j], T$i, max, nrow(L)),
               ns$tapply_full(wd$q_r$y[T$j], T$i, min, nrow(L)), ns$tapply_full(wd$q_r$y[T$j], T$i, max, nrow(L)))
  expect_equal(unname(bb), unname(ref))
  expect_true(all(is.na(bb[Matrix::rowSums(L > 0) == 0, ])))
  rows <- c(1, 5, 50, nrow(L))
  D <- Matrix::Diagonal(nrow(L), x = as.numeric(!(seq_len(nrow(L)) %in% rows)))
  expect_identical(ns$exclude_rows(L, rows), methods::as(Matrix::drop0(D %*% L), "dgCMatrix"))
  expect_identical(ns$exclude_rows(L, integer(0)), L)
})

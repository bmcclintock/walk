# Internal helpers for the `clipRast` argument of proc_Langevin() / proc_telem()

# Resolve `clipRast` to a window of whole raster cells (1-based rows/cols of `cell_data`).
#   clipRast = FALSE      -> NULL (no clipping)
#   clipRast = TRUE       -> bounding box of the error neighbourhoods (the `t_buf` region),
#                            i.e. [x - max_err, x + max_err] x [y - max_err, y + max_err] over all
#                            locations (locations with max_err <= 0 or NA contribute their point)
#   clipRast = number b   -> that bounding box expanded by b map units on every side
#                            (TRUE is the same as 0)
#   clipRast = SpatExtent -> that extent (anything terra::ext() accepts, e.g. a SpatRaster or
#                            SpatVector, also works)
# The window is snapped outward to the raster grid and intersected with the raster extent, so
# it equals the extent of terra::crop(cell_data, <window>, snap = "out").
#' @noRd
clip_window <- function(clipRast, cell_data, x, y, max_err){
  if(is.null(clipRast) || isFALSE(clipRast)) return(NULL)
  msg <- "'clipRast' must be FALSE, TRUE, a single non-negative buffer distance, or a terra::ext() object."
  data_box <- isTRUE(clipRast) || (is.numeric(clipRast) && is.null(dim(clipRast)) && !inherits(clipRast, "SpatExtent"))
  if(data_box){
    buf <- if(isTRUE(clipRast)) 0 else clipRast
    if(length(buf) != 1 || !is.finite(buf) || buf < 0) stop(msg, call. = FALSE)
    me <- ifelse(is.na(max_err) | max_err < 0, 0, max_err)
    ok <- is.finite(x) & is.finite(y)
    if(!any(ok)) stop("'clipRast': no finite locations to define the clipping region.", call. = FALSE)
    w <- c(xmin = min(x[ok] - me[ok]) - buf, xmax = max(x[ok] + me[ok]) + buf,
           ymin = min(y[ok] - me[ok]) - buf, ymax = max(y[ok] + me[ok]) + buf)
  } else if(is.logical(clipRast)){
    stop(msg, call. = FALSE)
  } else {
    e <- tryCatch(terra::ext(clipRast), error = function(err) NULL)
    if(is.null(e)) stop(msg, call. = FALSE)
    w <- as.vector(e)[c("xmin", "xmax", "ymin", "ymax")]
    if(any(!is.finite(w)) || w[["xmax"]] <= w[["xmin"]] || w[["ymax"]] <= w[["ymin"]])
      stop("'clipRast' extent is empty or not finite.", call. = FALSE)
  }
  re <- as.vector(terra::ext(cell_data))
  rs <- terra::res(cell_data); nr <- terra::nrow(cell_data); nc <- terra::ncol(cell_data)
  # column/row index containing a coordinate (clamped in double so huge values can't overflow)
  idx <- function(v, hi) as.integer(pmin(pmax(v, -1), hi + 1))
  if(data_box){
    # same closed-interval rule as the C++ neighbourhood search, so every cell that can touch
    # an error disc (even only at its edge) lies inside the window
    c0 <- floor((w[["xmin"]] - re[["xmin"]]) / rs[1]) + 1
    c1 <- floor((w[["xmax"]] - re[["xmin"]]) / rs[1]) + 1
    r0 <- floor((re[["ymax"]] - w[["ymax"]]) / rs[2]) + 1
    r1 <- floor((re[["ymax"]] - w[["ymin"]]) / rs[2]) + 1
  } else {
    # outward snap like terra::crop(snap = "out"); the tolerance stops an extent lying exactly
    # on cell edges from picking up an extra row/column
    tol <- 1e-8
    c0 <- floor((w[["xmin"]] - re[["xmin"]]) / rs[1] + tol) + 1
    c1 <- ceiling((w[["xmax"]] - re[["xmin"]]) / rs[1] - tol)
    r0 <- floor((re[["ymax"]] - w[["ymax"]]) / rs[2] + tol) + 1
    r1 <- ceiling((re[["ymax"]] - w[["ymin"]]) / rs[2] - tol)
  }
  c0 <- max(1L, idx(c0, nc)); c1 <- min(nc, idx(c1, nc))
  r0 <- max(1L, idx(r0, nr)); r1 <- min(nr, idx(r1, nr))
  if(c1 < c0 || r1 < r0) stop("'clipRast' extent does not overlap 'cell_data'.", call. = FALSE)
  list(r0 = r0, r1 = r1, c0 = c0, c1 = c1, nc_full = nc, nr_full = nr, data_box = data_box,
       nrow = r1 - r0 + 1L, ncol = c1 - c0 + 1L,
       ext = terra::ext(re[["xmin"]] + (c0 - 1) * rs[1], re[["xmin"]] + c1 * rs[1],
                        re[["ymax"]] - r1 * rs[2],       re[["ymax"]] - (r0 - 1) * rs[2]))
}

# Cell numbers (in the full `cell_data`) of every cell in the window, in clipped-raster order
#' @noRd
clip_cell_index <- function(win){
  rows <- win$r0:win$r1; cols <- win$c0:win$c1
  as.vector(t(outer((rows - 1) * win$nc_full, cols, `+`)))
}

# Treat cells outside the window as inadmissible (like NA cells), so normalisation, truncation
# and the "no admissible mass" checks work exactly as for masked cells.
#' @noRd
clip_mask <- function(mask, win){
  idx <- clip_cell_index(win)
  m <- rep(TRUE, length(mask))
  m[idx] <- mask[idx]
  m
}

# Convert full-raster cell numbers to clipped-raster cell numbers
#' @noRd
clip_renumber <- function(cell, win){
  row <- (cell - 1) %/% win$nc_full + 1
  col <- (cell - 1) %%  win$nc_full + 1
  (row - win$r0) * win$ncol + (col - win$c0 + 1)
}

# Number of locations whose error neighbourhood is cut by the window (sides of the window that
# coincide with the raster edge are ignored, since nothing exists beyond them anyway)
#' @noRd
clip_n_cut <- function(win, cell_data, x, y, max_err){
  e <- as.vector(win$ext); re <- as.vector(terra::ext(cell_data))
  lo_x <- if(win$c0 > 1) e[["xmin"]] else -Inf
  hi_x <- if(win$c1 < win$nc_full) e[["xmax"]] else Inf
  lo_y <- if(win$r1 < terra::nrow(cell_data)) e[["ymin"]] else -Inf
  hi_y <- if(win$r0 > 1) e[["ymax"]] else Inf
  me <- ifelse(is.na(max_err) | max_err < 0, 0, max_err)
  sum(x - me < lo_x | x + me > hi_x | y - me < lo_y | y + me > hi_y, na.rm = TRUE)
}

# make_walk_data(): bring 'cell_data' (and 'rast_mask') onto the clipped grid of proc_data$L.
# Accepts either the raster the data were processed on (it is cropped to clip$ext) or a raster
# that is already that grid (e.g. terra::crop(cell_data, terra::ext(clip$ext))). 'rast_mask' may
# likewise be on the full or the clipped grid.
#' @noRd
clip_apply <- function(clip, cell_data, rast_mask=NULL){
  ce <- terra::ext(clip$ext)
  on_grid <- function(r){
    tol <- 1e-6 * min(terra::res(r))
    all(dim(r)[1:2] == clip$dim) && all(abs(as.vector(terra::ext(r)) - as.vector(ce)) < tol)
  }
  to_grid <- function(r, what){
    if(on_grid(r)) return(r)
    # full raster: take the window's rows/columns exactly (no floating-point snapping)
    rc <- NULL
    if(terra::ncell(r) >= max(clip$cell)){
      rcl <- terra::rowColFromCell(r, clip$cell[c(1, length(clip$cell))])
      if(rcl[2,1] - rcl[1,1] + 1 == clip$dim[[1]] && rcl[2,2] - rcl[1,2] + 1 == clip$dim[[2]])
        rc <- r[rcl[1,1]:rcl[2,1], rcl[1,2]:rcl[2,2], drop=FALSE]
    }
    if(is.null(rc) || !on_grid(rc))
      stop("'", what, "' is neither the raster passed to proc_telem()/proc_Langevin() nor its ",
           "'clipRast' window (", paste(clip$dim, collapse=" x "), " cells with extent ",
           "'proc_data$clip$ext').", call.=FALSE)
    rc
  }
  cell_data <- to_grid(cell_data, "cell_data")
  if(!is.null(rast_mask)) rast_mask <- to_grid(rast_mask, "rast_mask")
  list(cell_data=cell_data, rast_mask=rast_mask)
}

# (one constructor, so clip elements made here and by clip_union() are identical())
#' @noRd
clip_make <- function(ext, win){
  list(ext = ext, dim = c(nrow = as.integer(win$nrow), ncol = as.integer(win$ncol)),
       cell = as.numeric(clip_cell_index(win)),
       full_dim = c(nrow = as.integer(win$nr_full), ncol = as.integer(win$nc_full)))
}

# The `clip` element stored in proc_* output
#' @noRd
clip_info <- function(win) clip_make(as.vector(win$ext), win)

# bind_proc_data(): put inputs clipped to different windows of the same raster onto one grid.
# `cl` is a list of `clip` elements (NULL = unclipped, i.e. the whole raster) and `nc` the number
# of columns of each input's L. Returns NULL if all windows are identical (nothing to do), else
# list(clip = <clip of the union window, or NULL for the whole raster>,
#      ncell = number of columns of the combined L,
#      map  = list of integer vectors: new column of each old column, per input).
#' @noRd
clip_union <- function(cl, nc){
  if(all(vapply(cl, function(z) identical(z, cl[[1]]), TRUE))) return(NULL)
  bad <- "All inputs must have been processed with the same 'cell_data'"
  hasc <- !vapply(cl, is.null, TRUE)
  fd <- lapply(cl[hasc], function(z) z$full_dim)
  if(any(vapply(fd, is.null, TRUE)))
    stop("Invalid 'clip' element (no 'full_dim'); use the output of proc_telem()/proc_Langevin().", call. = FALSE)
  if(!all(vapply(fd, function(z) identical(as.integer(z), as.integer(fd[[1]])), TRUE)))
    stop(bad, " (the clipped inputs come from rasters of different dimensions).", call. = FALSE)
  fd <- as.integer(fd[[1]])
  ntot <- fd[1] * fd[2]
  if(any(nc[!hasc] != ntot)) stop(bad, " (different numbers of cells).", call. = FALSE)
  # extents must lie on one grid: same resolution and the same origin
  rs <- t(vapply(cl[hasc], function(z) c((z$ext[2] - z$ext[1]) / z$dim[[2]], (z$ext[4] - z$ext[3]) / z$dim[[1]]), c(0, 0)))
  rc <- lapply(cl[hasc], function(z){ r <- (z$cell[c(1, length(z$cell))] - 1) %/% fd[2] + 1
                                       c <- (z$cell[c(1, length(z$cell))] - 1) %%  fd[2] + 1
                                       c(r0 = r[1], r1 = r[2], c0 = c[1], c1 = c[2]) })
  org <- t(vapply(seq_along(rc), function(k){ z <- cl[hasc][[k]]
    c(z$ext[1] - (rc[[k]][["c0"]] - 1) * rs[k, 1], z$ext[4] + (rc[[k]][["r0"]] - 1) * rs[k, 2]) }, c(0, 0)))
  tol <- 1e-6 * min(rs)
  if(any(abs(sweep(rs, 2, rs[1, ])) > tol) || any(abs(sweep(org, 2, org[1, ])) > tol))
    stop(bad, " (the clipped inputs are not on the same raster grid).", call. = FALSE)
  if(!all(hasc)){
    # an unclipped input: combine on the whole raster
    return(list(clip = NULL, ncell = ntot, map = lapply(seq_along(cl), function(k)
      if(is.null(cl[[k]])) seq_len(ntot) else as.integer(cl[[k]]$cell))))
  }
  r0s <- vapply(rc, `[[`, 1, "r0"); r1s <- vapply(rc, `[[`, 1, "r1")
  c0s <- vapply(rc, `[[`, 1, "c0"); c1s <- vapply(rc, `[[`, 1, "c1")
  r0 <- as.integer(min(r0s)); r1 <- as.integer(max(r1s))
  c0 <- as.integer(min(c0s)); c1 <- as.integer(max(c1s))
  win <- list(r0 = r0, r1 = r1, c0 = c0, c1 = c1, nc_full = fd[2], nr_full = fd[1],
              nrow = r1 - r0 + 1L, ncol = c1 - c0 + 1L)
  # take each side from an input window that has it, so the extent is exactly the one
  # proc_telem()/proc_Langevin() would give for this window
  ec <- lapply(cl[hasc], function(z) unname(z$ext))
  ext <- c(xmin = ec[[which.min(c0s)]][1], xmax = ec[[which.max(c1s)]][2],
           ymin = ec[[which.max(r1s)]][3], ymax = ec[[which.min(r0s)]][4])
  out <- clip_make(ext, win)
  list(clip = out, ncell = length(out$cell), map = lapply(cl, function(z) as.integer(clip_renumber(z$cell, win))))
}

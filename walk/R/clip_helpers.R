# Internal helpers for the `clipRast` argument of proc_Langevin() / proc_telem()

# Resolve `clipRast` to a window of whole raster cells (1-based rows/cols of `cell_data`).
#   clipRast = FALSE      -> NULL (no clipping)
#   clipRast = TRUE       -> bounding box of the error neighbourhoods (the `t_buf` region),
#                            i.e. [x - max_err, x + max_err] x [y - max_err, y + max_err] over all
#                            locations (locations with max_err <= 0 or NA contribute their point)
#   clipRast = SpatExtent -> that extent (anything terra::ext() accepts, e.g. a SpatRaster or
#                            SpatVector, also works)
# The window is snapped outward to the raster grid and intersected with the raster extent, so
# it equals the extent of terra::crop(cell_data, <window>, snap = "out").
#' @noRd
clip_window <- function(clipRast, cell_data, x, y, max_err){
  if(is.null(clipRast) || isFALSE(clipRast)) return(NULL)
  if(isTRUE(clipRast)){
    me <- ifelse(is.na(max_err) | max_err < 0, 0, max_err)
    ok <- is.finite(x) & is.finite(y)
    if(!any(ok)) stop("'clipRast = TRUE': no finite locations to define the clipping region.", call. = FALSE)
    w <- c(xmin = min(x[ok] - me[ok]), xmax = max(x[ok] + me[ok]),
           ymin = min(y[ok] - me[ok]), ymax = max(y[ok] + me[ok]))
  } else if(is.logical(clipRast)){
    stop("'clipRast' must be TRUE, FALSE, or a terra::ext() object.", call. = FALSE)
  } else {
    e <- tryCatch(terra::ext(clipRast), error = function(err) NULL)
    if(is.null(e)) stop("'clipRast' must be TRUE, FALSE, or a terra::ext() object.", call. = FALSE)
    w <- as.vector(e)[c("xmin", "xmax", "ymin", "ymax")]
    if(any(!is.finite(w)) || w[["xmax"]] <= w[["xmin"]] || w[["ymax"]] <= w[["ymin"]])
      stop("'clipRast' extent is empty or not finite.", call. = FALSE)
  }
  re <- as.vector(terra::ext(cell_data))
  rs <- terra::res(cell_data); nr <- terra::nrow(cell_data); nc <- terra::ncol(cell_data)
  # column/row index containing a coordinate (clamped in double so huge values can't overflow)
  idx <- function(v, hi) as.integer(pmin(pmax(v, -1), hi + 1))
  if(isTRUE(clipRast)){
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
  list(r0 = r0, r1 = r1, c0 = c0, c1 = c1, nc_full = nc,
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

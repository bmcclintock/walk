#' @title Map walk results (utilisation distribution or location predictions) to a raster
#' @description Turns the output of \code{\link{get_lim_ud}} or \code{\link{predict_ctmc}} into a
#' `SpatRaster` on the grid of the `cell_data` raster given to \code{\link{make_walk_data}}
#' (the full raster, also when `clipRast` was used). Cells that are not
#' part of the model (masked, `NA`, outside the window) are `NA`.
#' @param x A data frame from \code{\link{get_lim_ud}} (columns `cell`, `ud` and any `hpd_*`),
#' or a list from \code{\link{predict_ctmc}}.
#' @param rows For \code{\link{predict_ctmc}} output: which rows of `x$times` /
#' `x$local_state_prob` to map (one layer each). Defaults to all rows, which can be many layers.
#' @param template Optional `SpatRaster` with the grid of the `cell_data` given to
#' \code{\link{make_walk_data}}. Only needed if `x` does not carry the grid (the attribute is
#' dropped by some operations, e.g. `merge()`).
#' @return A `SpatRaster`. For a UD: layers `ud` and the `hpd_*` columns. For predictions:
#' one layer per row, named by the time stamp (and `id` for multi-individual data).
#' @examples
#' \dontrun{
#' ud <- get_lim_ud(fit, hpd = 0.95)
#' terra::plot(walk_rast(ud))
#' pr <- predict_ctmc(fit)
#' terra::plot(walk_rast(pr, rows = 1:4))
#' }
#' @export
walk_rast <- function(x, rows = NULL, template = NULL){
  is_ud <- is.data.frame(x) && all(c("cell", "ud") %in% names(x))
  is_pr <- !is.data.frame(x) && is.list(x) && !is.null(x$local_state_prob) && !is.null(x$cell_map)
  if(!is_ud && !is_pr) stop("'x' must be the output of get_lim_ud() or predict_ctmc().")
  if(is.null(template)) template <- attr(x, "cell_data")
  if(is.null(template)) stop("'x' does not carry the raster grid; supply 'template' (the 'cell_data' raster given to make_walk_data()).")
  if(inherits(template, "PackedSpatRaster")) template <- terra::unwrap(template)
  if(!inherits(template, "SpatRaster")) stop("'template' must be a SpatRaster.")
  r0 <- terra::rast(template[[1]]); terra::values(r0) <- NA_real_
  if(is_ud){
    if(any(x$cell > terra::ncell(r0))) stop("Cell numbers in 'x' exceed the number of cells of 'template'.")
    v <- c("ud", grep("^hpd_", names(x), value = TRUE))
    out <- terra::rast(lapply(v, function(z){ r <- r0; r[x$cell] <- x[[z]]; r }))
    names(out) <- v
    return(out)
  }
  if(is_pr){
    P <- x$local_state_prob
    if(is.null(rows)) rows <- seq_len(nrow(P))
    if(is.logical(rows)) rows <- which(rows)
    if(any(rows < 1 | rows > nrow(P))) stop("'rows' must be row numbers of 'x$local_state_prob' (1..", nrow(P), ").")
    cell <- x$cell_map$cell[order(x$cell_map$cellx)]
    if(any(cell > terra::ncell(r0))) stop("Cell numbers in 'x' exceed the number of cells of 'template'.")
    M <- matrix(NA_real_, terra::ncell(r0), length(rows))
    M[cell, ] <- as.matrix(Matrix::t(P[rows, , drop = FALSE]))
    out <- terra::rast(r0, nlyrs = length(rows)); terra::values(out) <- M
    tm <- x$times
    nm <- as.character(tm$timestamp[rows])
    if(length(unique(tm$id)) > 1) nm <- paste(tm$id[rows], nm)
    names(out) <- make.unique(nm)
    return(out)
  }
}

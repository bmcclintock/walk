#' @title Produce design data for use in fitting MMPP movement models
#' @param proc_data A sparse matrix with rows corresponding to time steps and columns corresponding to cells. 
#' The entries are probabilities that the animal is located in the corresponding cell and time.
#' @param cell_data A `SpatRaster` object from the `terra` package. If `proc_data` was made with
#' `clipRast` (it has a `clip` element), this can be the raster given to `proc_telem()` /
#' `proc_Langevin()` or the clipped window of it; the full raster is cropped to the window.
#' The model domain (the cells of the CTMC, and so of the utilization distribution) is the
#' `clipRast` window, or all of `cell_data` without `clipRast`; to restrict it, use `clipRast`
#' in `proc_telem()`/`proc_Langevin()`.
#' Memory use here grows with the number of cells in the domain and is usually the binding
#' constraint for large rasters (much more than the size of `proc_data$L`).
#' @param grad A character vector of names of `cell_data` layers for which gradient covariates will be constructed.
#' @param rast_mask Raster mask for inaccessible cells when `cell_data` is of type `SpatRaster` from the `terra` package. This is ignored
#' if \code{cell_data} is an `POLYGON` data frame from the `sf` package. With clipped `proc_data`
#' it may be on the full or the clipped grid, like `cell_data`.
#' @param directions Neighborhood structure, one of `"rook"` or `"queen"`.
#' @param mask_zero Logical. If `TRUE` (default) cells where `rast_mask` is `0` are dropped as
#' well as cells where it is `NA`. A mask built with
#' `terra::rasterize(land, r, field = 0, background = 1)` contains no `NA`s, so with
#' `mask_zero = FALSE` it would drop nothing.
#' @param add_diff Logical (default `TRUE`). Add difference covariates `d_<layer>` (value in
#' the destination cell minus value in the origin cell) to `q_m` for every numeric layer of
#' `cell_data`, so movement formulas such as `~ d_depth` can be used directly. Existing columns
#' (e.g. a layer that is itself called `d_depth`) are not overwritten. Other covariates can be
#' added to `q_r`/`q_m` by hand afterwards.
#' @param exclude Optional observations to exclude: rows of `proc_data$times` (= its column
#' `obs`) or a logical vector with one value per row. They are labelled `"user"` in
#' `times$excluded`. See \code{\link{exclude_obs}}.
#' @param vmax Optional maximum plausible speed, in map units per `time_unit` of the data (e.g.
#' m/h). Observations that would need a faster movement from their neighbours (allowing for
#' location error) are excluded and labelled `"speed"`. See \code{\link{exclude_obs}} for the
#' filter.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @details Cells that are `NA` in the first layer of `cell_data` are inadmissible (as in
#' `proc_telem()`/`proc_Langevin()`) and are dropped like masked cells. In the output,
#' `q_r$cell`, `q_m$cell` and `q_m$from_cell` are cell numbers of the `cell_data` raster supplied
#' here (also stored, first layer only, in `cell_data`), even when `proc_data` was clipped with
#' `clipRast`, so `get_lim_ud()` and `predict_ctmc()` results can be
#' mapped with e.g. `u <- terra::rast(cell_data); u[ud$cell] <- ud$ud`. `cellx` is the index of
#' the cell among the retained cells (the columns of `L`).
#'
#' **Excluded observations.** `times$excluded` gives, for every observation, why it is excluded
#' (`NA` if it is not): `"no_cells"` if its likelihood is zero in every cell of the model domain
#' (e.g. its whole error neighbourhood is masked; these are excluded automatically, with a
#' warning), `"user"` for `exclude`, and `"speed"` for the `vmax` filter, applied in that order.
#' An excluded observation keeps its time but its row of `L` is empty. In
#' \code{\link{fit_ctmc}} this is exactly equivalent to deleting it and measuring the next time
#' step from the previous retained observation; \code{\link{predict_ctmc}} still predicts its
#' location from the movement model. All exclusions work the same way, whatever their
#' reason. Use \code{\link{exclude_obs}} to exclude more observations later (e.g. `fit$zero_obs`).
#' @import dplyr 
#' @importFrom Matrix rowSums
#' @importFrom igraph graph_from_data_frame
#' @export
make_walk_data <- function(proc_data, cell_data, grad=NULL, rast_mask=NULL, 
                           directions="rook", mask_zero=TRUE,
                           add_diff=TRUE, exclude=NULL, vmax=NULL, debug=0){
  
  cells <- obs <- timestamp <- quad <- NULL
  
  if(debug==1) browser()
  if(!is.list(proc_data) || is.data.frame(proc_data) || is.null(proc_data$L) || is.null(proc_data$times))
    stop("'proc_data' must be the output of proc_telem()/proc_Langevin() (or bind_proc_data()) with ",
         "return_type = \"sparse\" (default) or \"dense\": a list with elements 'L' and 'times'.")
  if(is.matrix(proc_data$L)) proc_data$L <- as(Matrix::Matrix(proc_data$L, sparse=TRUE), "dgCMatrix")
  # check 'exclude' / 'vmax' before the (possibly expensive) rate-matrix data are built
  N <- nrow(proc_data$L)
  if(!is.null(exclude)){
    if(is.logical(exclude)){
      if(length(exclude) != N || anyNA(exclude)) stop("A logical 'exclude' must have one non-missing value per row of 'proc_data$times'.")
    } else if(!is.numeric(exclude) || any(is.na(exclude) | exclude < 1 | exclude > N | exclude != round(exclude)))
      stop("'exclude' must be row numbers of 'proc_data$times' (1..", N, ") or a logical vector.")
  }
  if(!is.null(vmax) && (!is.numeric(vmax) || length(vmax) != 1 || !is.finite(vmax) || vmax <= 0))
    stop("'vmax' must be a single positive number.")
  cell_name <- NULL
  quad_data <- NULL
  
  if(debug==2) browser()
  
  if(inherits(cell_data, "SpatRaster")){
    cell_data_in <- cell_data   # cell numbers in the output refer to this raster
    # proc_telem()/proc_Langevin() with 'clipRast': L is indexed on the clipped grid, so crop
    # cell_data (and rast_mask) to it when the full raster is supplied
    if(!is.null(proc_data$clip)){
      cr <- clip_apply(proc_data$clip, cell_data, rast_mask)
      cell_data <- cr$cell_data; rast_mask <- cr$rast_mask
    }
    if(ncol(proc_data$L) != terra::ncell(cell_data)) stop("'proc_data' was not built on 'cell_data' (different number of cells).")
    if(!is.null(rast_mask) && mask_zero){
      rast_mask <- terra::ifel(rast_mask == 0, NA, rast_mask)
    }
    if(!is.null(grad) && !all(grad %in% names(cell_data))) stop("'grad' layer(s) not found in 'cell_data': ", paste(setdiff(grad, names(cell_data)), collapse=", "))
    q_list <- make_q_data_rast(cell_data, grad=grad, rast_mask=rast_mask, directions=directions,
                               grad_rast=if(terra::ncell(cell_data) != terra::ncell(cell_data_in)) cell_data_in else NULL)
    attr(q_list$q_m, "directions") <- directions
    attr(q_list$q_m, "cov_class") <- "SpatRaster"
  } else if(inherits(cell_data, "sf")){
    stop("'cell_data' must be a 'SpatRaster' from the terra package.")
    # q_list <- make_q_data_sf(cell_data, cell_name)
    # attr(q_list$q_m, "directions") <- NULL
    # attr(q_list$q_m, "cov_class") <- "sf"
  }
  out <- list()
  out$L <- proc_data$L[,q_list$q_r$cell]
  out$times <- proc_data$times
  if(is.null(out$times$seg)) out$times$seg <- 1L
  if(is.null(out$times$id)) out$times$id <- "1"
  if(is.null(attr(out$times, "time_unit"))) warning("'proc_data$times' has no 'time_unit' attribute; 'dt' units are unknown.")
  zl <- which(Matrix::rowSums(out$L) == 0)
  if(length(zl)){
    zid <- paste0(out$times$id[zl], ":", out$times$obs[zl])
    warning(length(zl), " observation(s) have likelihood 0 in every cell of the model domain (after ",
            "'rast_mask' and NA 'cell_data' cells are removed) and are excluded (times$excluded == \"no_cells\"). ",
            "(id:obs) ", paste(utils::head(zid, 10), collapse=", "), if(length(zid) > 10) ", ...",
            ". If that is unexpected, check 'rast_mask' and the location errors.", call. = FALSE)
  }
  # 'clipRast' works on a window of 'cell_data'. Report cell numbers in the
  # raster supplied by the user so that e.g. r[ud$cell] <- ud$ud works on that raster.
  if(inherits(cell_data, "SpatRaster")){
    if(terra::ncell(cell_data) != terra::ncell(cell_data_in) ||
       any(abs(as.vector(terra::ext(cell_data)) - as.vector(terra::ext(cell_data_in))) > 1e-6*min(terra::res(cell_data)))){
      orig <- terra::cellFromXY(cell_data_in, terra::xyFromCell(cell_data, seq_len(terra::ncell(cell_data))))
      q_list$q_r$cell <- orig[q_list$q_r$cell]
      q_list$q_m$cell <- orig[q_list$q_m$cell]
      q_list$q_m$from_cell <- orig[q_list$q_m$from_cell]
    }
    out$cell_data <- terra::wrap(cell_data_in[[1]])
  }
  out$q_r <- q_list$q_r
  out$q_m <- q_list$q_m
  # difference covariates d_<layer> = value in the destination cell - value in the origin cell
  if(isTRUE(add_diff)){
    lyr <- names(cell_data)
    lyr <- lyr[vapply(lyr, function(z) is.numeric(out$q_m[[z]]) && is.numeric(out$q_m[[paste0("from_", z)]]), logical(1))]
    lyr <- lyr[!(paste0("d_", lyr) %in% names(out$q_m))]   # never overwrite an existing column
    for(z in lyr) out$q_m[[paste0("d_", z)]] <- out$q_m[[z]] - out$q_m[[paste0("from_", z)]]
  }
  class(out) <- c(class(out), "walk_ddl")
  # excluded observations: "no_cells" (empty L rows), then 'exclude' ("user"), then 'vmax' ("speed")
  out <- exclude_obs(out, obs = exclude, vmax = vmax, reason = "user")
  # Check for irreducubility
  g <- graph_from_data_frame(d = out$q_m[,c("from_cellx","cellx")], directed = TRUE)
  scc <- components(g, mode = "strong")
  if(scc$no>1){
    warning("There appears to be ", scc$no, " disconnected groups of cells in the masked habitat.",
            " This can cause fit failure becuase the implied CTMC is not irreducible. ",
            "In addition, the limiting distribution does not exist.")
  }
  if(debug==3) browser()
  return(out)
}





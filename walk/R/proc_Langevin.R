#' @title Convert `simLangevin` or `dataLangevin` data to a likelihood surface on a `SpatRaster`
#' @param data A `simLangevin` or `dataLangevin` object from the `langevinSSM` package, or a
#' **list** of them (e.g. one per individual).
#' @param cell_data A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param return_type Type of object returned. One if `"data.frame"`, `"sparse"` (sparse matrix), 
#' `"matrix_df"` (matrix form of `"data.frame"`), or `"dense"` (dense matrix).
#' @param max_err The maximum error in meters. If unspecified it will be set to
#' 4 times the maximum error standard deviation as determined by the UERE and HDOP
#' of the telemetry data.
#' @param trunc The smallest probability value that is considered to be > 0. Defaults to 1.0e-8.
#' Each observation's cell probabilities are normalised to sum to 1 and values below `trunc` are
#' then set to 0 without renormalising; this only adds a parameter-free constant to the
#' log-likelihood.
#' @param  return_inadmissible Logical. If there are locations outside the admissible areas of movement return the observations and raster cell numbers.
#' @param nthreads Number of OpenMP threads used for the cell probability calculations. Defaults to 1.
#' Has no effect if `walk` was built without OpenMP (a warning is given); check with `walk:::omp_info()`.
#' @param neighborhood How the error neighbourhood of each location is found. `"disc"` (default, 
#' fastest) uses every raster cell that touches an exact circle of radius `max_err`, computed 
#' arithmetically from the raster geometry. `"buffer"` uses `terra::buffer()` + 
#' `terra::cells(touches=TRUE)`. The polygonised 
#' buffer sits slightly inside the true circle, so it can miss a few edge cells; otherwise the 
#' results are identical. Probabilities are computed in C++ either way.
#' @param clipRast The model domain: which part of `cell_data` the movement model (and so the
#' utilization distribution) is defined on. `FALSE` (default) uses the whole raster. A number
#' `b >= 0` uses the bounding box of the error neighbourhoods (the circles of radius `max_err`
#' around all locations) expanded by `b` map units on every side; `TRUE` is the same as `0`
#' (no buffer). A `terra::ext()` object (or anything `terra::ext()` accepts, such as a `SpatRaster`
#' or `SpatVector`) uses that extent. The region is snapped outward to whole cells of
#' `cell_data` and limited to the raster, so it matches
#' `terra::crop(cell_data, <region>, snap = "out")`. Cells outside the region are treated like
#' `NA` cells: each location's likelihood is renormalised over the cells inside it. A warning is
#' given if a user-supplied extent cuts through any location's error neighbourhood. When clipping,
#' `L` has one column per cell of the clipped raster (cell numbers refer to the clipped raster)
#' and the output gains a `clip` element (an attribute for the data.frame/matrix return types)
#' with `ext` (the snapped extent as `c(xmin, xmax, ymin, ymax)`), `dim` (rows and columns of the
#' clipped raster), `cell` (cell number in `cell_data` of each column of `L`) and `full_dim`
#' (rows and columns of `cell_data`). `make_walk_data()` can be given the full `cell_data` (and a
#' full-size `rast_mask`): it crops them to the clipped grid itself, and its output reports cell
#' numbers of the full raster. `return_inadmissible` output always uses `cell_data` cell numbers.
#' Outputs with different regions (e.g. one call per individual with `clipRast = 5000`) can be
#' combined with \code{\link{bind_proc_data}}, which uses the bounding box of their regions.
#' See 'Model domain' below for how to choose the region.
#' @param id Name of the column in `data` identifying individuals (default `"id"`). If the column
#' is absent, list elements are named by their list names (or position), and a single object is
#' individual `"1"`. Time differences are computed **within** individual and every individual
#' starts a new, independent Markov chain in the likelihood. Rows of one individual must be
#' contiguous and in time order. Set to `NULL` to ignore any `id` column.
#' @param time_unit Unit for the stored time differences `times$dt`, and so for the fitted rates.
#' `NULL` (default) uses the time unit of `data`: the `"time.unit"` attribute that
#' `langevinSSM::formatData()` sets (and `simLangevin()` keeps when it simulates from data).
#' Otherwise one of `"secs"`, `"mins"`, `"hours"`, `"days"` or `"weeks"` (`"seconds"` and
#' `"minutes"` are also accepted), and `dt` is converted to that unit. If `data` has no
#' `"time.unit"` attribute: POSIXct/Date dates default to `"hours"` (the default `time.unit` of
#' `langevinSSM::formatData()`); numeric dates (e.g.
#' `simLangevin()` output simulated without data) are used as they are, in the units of the
#' `date` column, and the stored unit is `"native"` unless `time_unit` is given, in which case
#' it only names the unit of the numeric dates (no conversion).
#' @param max_gap Optional. Gaps longer than `max_gap` (in `time_unit`) start a new, independent
#' segment (the chain is restarted from `delta`). This removes the cost and the numerical
#' liability of propagating the CTMC across very long gaps without discarding any locations.
#' @details This function takes the error information in `data` to 
#' produce a likelihood surface over the `SpatRaster` defined by 
#' `cell_data` for each location. This can then be passed to 
#' CTMC HMM fitting functions. The arguments and the output have the same form as
#' \code{\link{proc_telem}}, so data processed by either function can be combined with
#' \code{\link{bind_proc_data}}. The `times` data frame returned with the likelihood
#' contains columns `obs`, `id`, `seg`, `timestamp` and `dt` (attributes `time_unit`, and
#' `ts_scale`, the factor converting differences of `as.numeric(timestamp)` to `dt`); `dt` is `0`
#' at the first observation of every segment. Observations with zero likelihood in every admissible
#' cell get an all-zero row in `L`; \code{\link{make_walk_data}} excludes them (`times$excluded == "no_cells"`).
#' 
#' @author Devin S. Johnson
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells res values
#' @importFrom methods as
#' @section Model domain:
#' The CTMC is defined on the cells of the domain only, and movement out of it is impossible
#' (the domain edge acts as a reflecting boundary). So the domain is part of the model, not just
#' a computational shortcut:
#' * Memory and run time of `make_walk_data()` and `fit_ctmc()` grow with the number of cells in
#'   the domain. For large rasters (millions of cells) `make_walk_data()` on the full raster can
#'   exhaust memory, so restrict the domain here with `clipRast`. This is the only place to do so:
#'   doing it here also keeps `L` small.
#' * The utilization distribution (`get_lim_ud()`, `walk_rast()`) is the stationary distribution
#'   on the domain. With `clipRast = TRUE` it ends at the edge of the location error
#'   neighbourhoods; use a buffer (`clipRast = <distance>`) large enough to cover the area where
#'   you want the UD. With `norm = TRUE` (the default in `ctmc_model()`) the UD near the edge
#'   depends a little on where the edge is, and with `delta = "stationary"` the parameter
#'   estimates can depend on the domain too; increase the buffer until the results of interest
#'   no longer change.
#' * The log-likelihood (and so AIC) depends on the domain, e.g. with `delta = "uniform"` it
#'   contains `log(number of cells)` once per track segment. Only compare models fitted on the
#'   same domain.
#' @export
proc_Langevin <- function(data, cell_data, return_type="sparse", max_err=NULL, trunc=1.0e-8, 
                          return_inadmissible=FALSE, nthreads=1, neighborhood=c("disc","buffer"),
                          clipRast=FALSE, id="id", time_unit=NULL, max_gap=NULL){
  neighborhood <- match.arg(neighborhood)
  proc_check_return_type(return_type)
  proc_check_threads(nthreads)
  # Check arguments
  dl <- proc_as_list(data, c("simLangevin", "dataLangevin"),
                     "a 'simLangevin' or 'dataLangevin' object from the langevinSSM package")
  if(!inherits(cell_data, "SpatRaster")) stop("'cell_data' must be a terra::SpatRaster object!")
  xl <- dl$x
  
  # Individual / segment structure
  idl <- proc_ids(xl, id, dl$is_list)
  if(any(vapply(xl, function(d) is.null(d$date), TRUE))) stop("'data' has no 'date' column.")
  dcl <- unique(vapply(xl, function(d) if(inherits(d$date, "POSIXt")) "POSIXct" else class(d$date)[1], ""))
  if(length(dcl) > 1) stop("The 'date' columns of the inputs have different classes (", paste(dcl, collapse=", "), ").")
  ts <- do.call(c, lapply(xl, function(d) d$date))
  tr <- proc_resolve_time(ts, lapply(xl, proc_obj_unit), time_unit, default="hours")
  times <- proc_times(unlist(idl, use.names=FALSE), ts, tr, max_gap)
  # Stack only the columns used below (list elements may carry different extra columns)
  cols <- c("x", "y", "smaj", "smin", "eor", "x.err", "y.err")
  d <- do.call(rbind, lapply(xl, function(z){
    z <- as.data.frame(z)
    for(cc in setdiff(cols, names(z))) z[[cc]] <- NA_real_
    z[, cols]
  }))
  n <- nrow(d)
  
  # Evaluate error ellipse covariance
  cov_df <- argos_diag_to_cov(d$smaj, d$smin, 180*d$eor/pi)
  
  # x and y standard deviation errors (assumes no correlation, i.e., rho_o = 0)
  sderrInd <- which(!is.na(d$x.err) & !is.na(d$y.err))
  if(length(sderrInd)){
    rho_o <- 0
    v1 <- d$x.err[sderrInd]^2
    v2 <- d$y.err[sderrInd]^2
    c12 <- d$x.err[sderrInd] * d$y.err[sderrInd] * rho_o
    cov_df$cov.x.x[sderrInd] <- v1
    cov_df$cov.y.y[sderrInd] <- v2
    cov_df$cov.x.y[sderrInd] <- c12
    cov_df$cov.check[sderrInd] <-  (v1*v2-c12^2) > 0
  }
  
  x <- as.numeric(d$x)
  y <- as.numeric(d$y)
  
  if(is.null(max_err)) {
    max_err <- 4 * sqrt(pmax(cov_df$cov.x.x, cov_df$cov.y.y))
  } else if(length(max_err) == 1) {
    max_err <- rep(max_err, n)
  }
  
  proc_lik(x, y, cov_df$cov.x.x, cov_df$cov.y.y, cov_df$cov.x.y, max_err, cell_data, times,
           return_type, trunc, return_inadmissible, nthreads, neighborhood, clipRast)
}

#' @title Combine processed telemetry from several individuals or data sources
#' @description Stack the output of several calls to \code{\link{proc_telem}} and/or
#' \code{\link{proc_Langevin}} (e.g. one per individual, possibly run in parallel, or
#' ctmm telemetry for some animals and langevinSSM data for others) into a single object
#' that \code{\link{make_walk_data}} and \code{\link{fit_ctmc}} will treat as a joint,
#' multi-individual data set.
#' @param ... Objects returned by `proc_telem(..., return_type="sparse")` or
#' `proc_Langevin(..., return_type="sparse")`, in any combination, or a single list of them.
#' @param ids Optional character vector of individual names, one per input. Overrides the
#' `id` column of each input's `times` (each input is then treated as one individual).
#' @details All inputs must have been built on the same `cell_data`. Inputs made with different
#' `clipRast` regions (e.g. `clipRast = TRUE` or a buffer distance per individual) are combined on
#' the smallest window of `cell_data` that contains all of them (the bounding box of the
#' individual windows, so the model domain is the same as if all individuals had been processed
#' in one call with that window); if any input is unclipped the result is unclipped. All inputs must use the same
#' `time_unit` and the same kind of time stamps (all POSIXct, or all numeric in the same units).
#' Segments are renumbered so they stay distinct, and an individual may not appear
#' in more than one input (use `max_gap`, not separate inputs, to split one animal's track).
#' @return A list with `L` (sparse), `times` (columns `obs`, `id`, `seg`, `timestamp`, `dt`,
#' attributes `time_unit` and `ts_scale`) and, if the inputs were clipped, `clip`; the same form as the output of
#' \code{\link{proc_telem}} and \code{\link{proc_Langevin}}.
#' @export
bind_proc_data <- function(..., ids=NULL){
  x <- list(...)
  if(length(x)==1 && is.list(x[[1]]) && is.null(x[[1]]$L)) x <- x[[1]]
  if(length(x)==0) stop("No inputs.")
  if(!all(vapply(x, function(z) is.list(z) && (inherits(z$L, "Matrix") || is.matrix(z$L)) && is.data.frame(z$times), TRUE)))
    stop("All inputs must be output of proc_telem() or proc_Langevin() with return_type = \"sparse\" or \"dense\".")
  if(!is.null(ids) && length(ids) != length(x)) stop("'ids' must have one entry per input.")
  nc <- vapply(x, function(z) ncol(z$L), 1)
  cl <- lapply(x, function(z) z$clip)
  # inputs clipped to different windows of the same raster are put on the union window
  un <- clip_union(cl, nc)
  if(!is.null(un)){
    for(k in seq_along(x)){
      if(is.null(cl[[k]])) next   # already on the whole raster
      Lk <- as(as(x[[k]]$L, "CsparseMatrix"), "generalMatrix")
      M <- Matrix::sparseMatrix(i = seq_len(ncol(Lk)), j = un$map[[k]], x = 1, dims = c(ncol(Lk), un$ncell))
      x[[k]]$L <- Lk %*% M
    }
    cl <- rep(list(un$clip), length(x))
    nc <- vapply(x, function(z) ncol(z$L), 1)
  }
  if(length(unique(nc))!=1) stop("All inputs must have been processed with the same 'cell_data' (different numbers of cells).")
  tu <- unique(vapply(x, function(z){ u <- attr(z$times, "time_unit"); if(is.null(u)) NA_character_ else u }, ""))
  if(length(tu)!=1 || is.na(tu)) stop("All inputs must have the same (non-missing) 'time_unit' attribute.")
  tsc <- vapply(x, function(z) proc_ts_scale(z$times), 1)
  if(any(abs(tsc - tsc[1]) > 1e-12 * abs(tsc[1])))
    stop("The inputs have time stamps on different scales (e.g. POSIXct and numeric dates). ",
         "Use the same kind of time stamps for every input.")
  seg_off <- 0L
  tl <- vector("list", length(x))
  for(k in seq_along(x)){
    tk <- x[[k]]$times
    if(nrow(tk) != nrow(x[[k]]$L)) stop("Input ", k, ": 'times' and 'L' have different numbers of rows.")
    if(is.null(tk$seg)) tk$seg <- 1L
    if(is.null(tk$id)) tk$id <- as.character(k)
    if(!is.null(ids)) tk$id <- as.character(ids[k])
    tk$id <- as.character(tk$id)
    tk$dt[1] <- 0
    tk$seg <- as.integer(factor(tk$seg, levels=unique(tk$seg))) + seg_off
    seg_off <- max(tk$seg)
    tl[[k]] <- tk[, c("obs", "id", "seg", "timestamp", "dt")]
  }
  times <- do.call(rbind, tl)
  rownames(times) <- NULL
  id_in <- unlist(lapply(tl, function(tk) unique(tk$id)), use.names=FALSE)
  if(any(duplicated(id_in))) stop("The same individual id appears in more than one input: ",
                                  paste(unique(id_in[duplicated(id_in)]), collapse=", "),
                                  ". Use 'ids' to rename, or process each animal in one call (use 'max_gap' to split a track).")
  times$obs <- seq_len(nrow(times))
  attr(times, "time_unit") <- tu
  attr(times, "ts_scale") <- tsc[1]
  L <- do.call(rbind, lapply(x, function(z) as(z$L, "CsparseMatrix")))
  out <- list(L=as(as(L, "generalMatrix"), "dgCMatrix"), times=times)
  if(!is.null(cl[[1]])) out$clip <- cl[[1]]
  out
}

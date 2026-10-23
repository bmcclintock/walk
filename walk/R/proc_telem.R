#' @title Convert `ctmm` telemetry data object to a `SpatRaster` stack from the `terra` package.
#' @param data A `telemetry` object from the `ctmm` package
#' @param cell_data A `\link[terra]{SpatRaster}` stack of covariates that will be used in CTMC 
#' movement modeling. Cells with `NA` values will be considered areas to the animal cannot travel, i.e., 
#' likelihood surfaces will be `0` for those cells.
#' @param time_unit Unit of measurement for numeric conversion of the POSIX timestamps 
#' in the telemetry data. Defaults to `"hours"`. Can specify, e.g., 
#'  `"seconds"` or `"days"`
#' @param return_type Type of object returned. One if `"data.frame"`, `"sparse"` (sparse matrix), 
#' `"matrix_df"` (matrix form of `"data.frame"`), or `"dense"` (dense matrix).
#' @param max_err The maximum error in meters. If unspecified it will be set to
#' 4 times the maximum error standard deviation as determined by the UERE and HDOP
#' of the telemetry data.
#' @param trunc The smallest probability value that is considered to be > 0. Defaults to 1.0e-8.
#' @param  return_inadmissible Logical. If there are locations outside the admissible areas of movement return the observations and raster cell numbers.
#' @param nthreads Number of OpenMP threads used for the cell probability calculations. Defaults to 1.
#' Has no effect if `walk` was built without OpenMP (a warning is given); check with `walk:::omp_info()`.
#' @param neighborhood How the error neighbourhood of each location is found. `"disc"` (default) 
#' uses exactly the raster cells that intersect the circle of radius `max_err`, computed 
#' arithmetically from the raster geometry. `"buffer"` uses `terra::buffer()` + 
#' `terra::cells(touches=TRUE)`; the buffer is a 40-sided polygon 
#' inscribed in the circle, so it can miss a few cells at the edge.
#' @param clipRast Restrict the returned likelihood matrix to part of `cell_data`. `FALSE` (default)
#' keeps the full raster, so `L` has `terra::ncell(cell_data)` columns. `TRUE` clips to the
#' error-buffer (`t_buf`) region: the bounding box of the circles of radius `max_err` around all
#' locations. A `terra::ext()` object (or anything `terra::ext()` accepts, such as a `SpatRaster`
#' or `SpatVector`) clips to that extent. The region is snapped outward to whole cells of
#' `cell_data` and limited to the raster, so it matches
#' `terra::crop(cell_data, <region>, snap = "out")`. Cells outside the region are treated like
#' `NA` cells: each location's likelihood is renormalised over the cells inside it. A warning is
#' given if a user-supplied extent cuts through any location's error neighbourhood. When clipping,
#' `L` has one column per cell of the clipped raster (cell numbers refer to the clipped raster)
#' and the output gains a `clip` element (an attribute for the data.frame/matrix return types)
#' with `ext` (the snapped extent as `c(xmin, xmax, ymin, ymax)`), `dim` (rows and columns of the
#' clipped raster) and `cell` (cell number in `cell_data` of each column of `L`). Use
#' `terra::crop(cell_data, terra::ext(out$clip$ext))` to get the matching covariate raster for
#' `make_walk_data()`. `return_inadmissible` output always uses `cell_data` cell numbers.
#' @details This function takes the HDOP information in the `telemetry` object to 
#' produce a `SpatRaster` likelihood surface over the `SpatRaster` defined by 
#' the `raster` argument for each location. This can then be passed to 
#' CTMC HMM fitting functions. 
#' 
#' @author Devin S. Johnson
#' @importFrom Matrix sparseMatrix
#' @importFrom terra vect buffer cells crds is.lonlat crs project res values
#' @importFrom ctmm uere
#' @importFrom methods as
#' @export
proc_telem <- function(data, cell_data, time_unit="hours", return_type="sparse", max_err=NULL, trunc=1.0e-8, 
                       return_inadmissible=FALSE, nthreads=1, neighborhood=c("disc","buffer"),
                       clipRast=FALSE){
  neighborhood <- match.arg(neighborhood)
  if(nthreads > 1){
    oi <- omp_info()
    if(!oi$enabled){
      warning("'walk' was compiled without OpenMP support, so 'nthreads' is ignored and 1 thread is used. ",
              "See walk:::omp_info().", call.=FALSE)
    } else if(oi$thread_limit < nthreads){
      warning("OpenMP thread limit (OMP_THREAD_LIMIT) is ", oi$thread_limit, "; using at most that many threads.", call.=FALSE)
    }
  }
  # Check arguments
  if(!inherits(data, "telemetry")) stop("'data' must be a ctmm::telemery object!")
  if(!inherits(cell_data, "SpatRaster")) stop("'cell_data' must be a terra::SpatRaster object!")
  if(is.lonlat(cell_data)) stop("'cell_data' must be projected with units in meters.")
  
  n <- nrow(data)
  
  # Evaluate error covariance
  err_nms <- c("COV.x.x", "COV.x.y", "COV.y.y")
  if(all(err_nms %in% colnames(data))){
    cov_xx <- data$COV.x.x
    cov_yy <- data$COV.y.y
    cov_xy <- data$COV.x.y
  } else {
    if(! "HDOP" %in% colnames(data)){
      warning("'HDOP' not detected in data it will be set to 1.")
      data$HDOP <- 1
    }
    u <- summary(uere(data))[,,"horizontal"]
    cov_xx <- (data$HDOP*u['est'])^2/2
    cov_yy <- (data$HDOP*u['est'])^2/2
    cov_xy <- rep(0, length(cov_xx))
  }
  cov_xx <- as.numeric(cov_xx); cov_yy <- as.numeric(cov_yy); cov_xy <- as.numeric(cov_xy)
  
  # Project locations to the raster CRS
  cell_data_prj <- crs(cell_data, proj=TRUE)
  telem_pts <- terra::vect(cbind(data$longitude, data$latitude), crs="epsg:4326") |>
    terra::project(cell_data_prj)
  xy <- terra::crds(telem_pts)
  
  if(is.null(max_err)) max_err <- 4*sqrt(pmax(cov_xx,cov_yy))
  max_err <- rep_len(as.numeric(max_err), n)
  if(any(is.na(max_err) | max_err <= 0)) stop("'max_err' must be positive (check the error covariance / HDOP values).")
  
  # Raster geometry and inadmissible-cell mask (NA in first layer, as in the original)
  e <- as.vector(terra::ext(cell_data))
  geom <- c(e[["xmin"]], e[["ymax"]], terra::res(cell_data), 
            terra::nrow(cell_data), terra::ncol(cell_data))
  mask <- is.na(terra::values(cell_data[[1]], mat=FALSE))
  
  # Optional clipping region: cells outside it are treated as inadmissible (like NA cells)
  win <- clip_window(clipRast, cell_data, xy[,1], xy[,2], max_err)
  if(!is.null(win)){
    mask <- clip_mask(mask, win)
    if(!isTRUE(clipRast)){
      n_cut <- clip_n_cut(win, cell_data, xy[,1], xy[,2], max_err)
      if(n_cut > 0) warning(n_cut, " location(s) have error neighbourhoods extending beyond 'clipRast'; ",
                            "likelihood mass outside it is dropped and the rest renormalised.", call.=FALSE)
    }
  }

  # Cell likelihoods: computed, masked, normalised and truncated in C++
  if(neighborhood == "disc"){
    res <- bvn_lik_disc(xy[,1], xy[,2], cov_xx, cov_yy, cov_xy, max_err, geom, mask, 
                        trunc, as.integer(nthreads))
  } else {
    t_buf <- terra::buffer(telem_pts, max_err)
    extr <- terra::cells(cell_data, t_buf, touches=TRUE)
    o <- order(extr[,'ID'])
    res <- bvn_lik_list(as.integer(extr[o,'ID']), as.numeric(extr[o,'cell']), xy[,1], xy[,2], 
                        cov_xx, cov_yy, cov_xy, max_err, geom, mask, trunc, as.integer(nthreads))
  }
  
  zero_ind <- length(res$inad_obs) > 0
  if(return_inadmissible){
    if(!zero_ind) return(NULL)
    return(cbind(obs=res$inad_obs, cell=res$inad_cell, lik=res$inad_lik))
  }
  if(zero_ind) warning("Some observations have likelihood values of 0 in all 'cell_data' cells!")
  
  times <- data.frame(obs=1:n, timestamp = data$timestamp)
  dt <- diff(times$timestamp) |> `units<-`(time_unit)
  times <- cbind(times, dt=c(0,dt))
  attr(times, "time_unit") <- time_unit
  
  ncells <- prod(dim(cell_data)[1:2])
  clip <- NULL
  if(!is.null(win)){
    res$cell <- clip_renumber(res$cell, win)
    ncells <- win$nrow * win$ncol
    clip <- list(ext = as.vector(win$ext), dim = c(nrow = win$nrow, ncol = win$ncol),
                 cell = clip_cell_index(win))
  }

  if(return_type %in% c("sparse", "dense")){
    M <- Matrix::sparseMatrix(i = res$obs, j = res$cell, x = res$lik, dims = c(n, ncells))
    if(return_type == "sparse") M <- as(M, "dgCMatrix") else M <- as.matrix(M)
    out <- list(L=M, times=times)
    if(!is.null(clip)) out$clip <- clip
    return(out)
  }
  out <- cbind(obs=res$obs, cell=res$cell, lik=res$lik)
  rm(res)
  if(!is.null(clip)) attr(out, "clip") <- clip
  if(return_type=="data.frame"){
    out <- as.data.frame(out)
    out <- merge(out, times, by="obs")
    if(!is.null(clip)) attr(out, "clip") <- clip
    return(out)
  } else if(return_type=="matrix_df"){
    return(out)
  } else{
    warning("Unknown 'return_type' returning 'matrix_df'")
    return(out)
  }
}

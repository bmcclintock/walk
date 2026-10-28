# Internal helpers shared by proc_telem() and proc_Langevin().
# Both functions only differ in how they read locations / error covariances
# from their input class; everything else (individual/segment bookkeeping,
# cell likelihoods, clipping and output assembly) lives here so the two
# behave identically.

# Seconds per 'time_unit' (difftime unit names, plus long-form aliases).
proc_time_unit <- function(time_unit){
  if(!is.character(time_unit) || length(time_unit) != 1) stop("'time_unit' must be a single character string.")
  tu <- switch(time_unit, secs=, seconds="secs", mins=, minutes="mins", hours="hours",
               days="days", weeks="weeks",
               stop("'time_unit' must be one of 'secs', 'mins', 'hours', 'days', 'weeks'."))
  list(unit=tu, mult=c(secs=1, mins=60, hours=3600, days=86400, weeks=604800)[[tu]])
}

# Time unit stored on an input object (langevinSSM: attr(, "time.unit")), normalised, or NULL.
proc_obj_unit <- function(d){
  u <- attr(d, "time.unit", exact=TRUE)
  if(is.null(u)) u <- attr(d, "time_unit", exact=TRUE)
  if(is.null(u) || length(u) != 1 || is.na(u)) return(NULL)
  tryCatch(proc_time_unit(as.character(u))$unit, error=function(e)
    stop("Unrecognised time unit attribute '", u, "' on 'data'. Set 'time_unit' explicitly.", call.=FALSE))
}

# Resolve the time unit of the output and the factor that converts differences of
# as.numeric(timestamp) into that unit.
#   ts:        stacked time stamps (POSIXct, Date or numeric)
#   obj_units: list of per-input unit attributes (NULL entries allowed)
#   time_unit: user value (NULL = use the object's unit)
#   default:   unit used for POSIXct/Date time stamps when neither time_unit nor an object
#              unit is given ("hours", the default time.unit of langevinSSM::formatData())
# Rules:
#   * POSIXct/Date time stamps: output unit = time_unit, else the object's unit, else 'default'.
#   * numeric time stamps with an object unit: they are in that unit; converted to time_unit if given.
#   * numeric time stamps without an object unit (e.g. simLangevin output): they are used as they
#     are; time_unit, if given, only names their unit; otherwise the unit is "native"
#     (the units of the 'date' column).
proc_resolve_time <- function(ts, obj_units, time_unit, default){
  ou <- unique(unlist(lapply(obj_units, function(u) if(is.null(u)) NA_character_ else u)))
  if(length(ou) > 1) stop("The inputs have different time units (", paste(ou, collapse=", "),
                          "). Convert them to a common unit first.", call.=FALSE)
  obj <- if(length(ou) == 0 || is.na(ou)) NULL else ou
  if(inherits(ts, "POSIXt")) ts_secs <- 1
  else if(inherits(ts, "Date")) ts_secs <- 86400
  else if(is.numeric(ts)) ts_secs <- if(is.null(obj)) NA_real_ else proc_time_unit(obj)$mult
  else stop("Time stamps must be POSIXct, Date or numeric.", call.=FALSE)
  unit <- if(!is.null(time_unit)) proc_time_unit(time_unit)$unit else if(!is.null(obj)) obj else
    if(is.na(ts_secs)) "native" else proc_time_unit(default)$unit
  scale <- if(unit == "native" || is.na(ts_secs)) 1 else ts_secs / proc_time_unit(unit)$mult
  list(unit=unit, scale=scale)
}

# Factor converting differences of as.numeric(times$timestamp) into times$dt units.
proc_ts_scale <- function(times){
  s <- attr(times, "ts_scale", exact=TRUE)
  if(is.null(s)) stop("'times' has no 'ts_scale' attribute. Re-run proc_telem() or proc_Langevin().", call.=FALSE)
  s
}

# Split a list of objects (one per individual) from a single object.
# Returns a list of objects; 'is_list' tells whether the input was a list.
proc_as_list <- function(data, cls, what){
  if(inherits(data, cls)) return(list(x=list(data), is_list=FALSE))
  if(is.list(data) && !is.data.frame(data) && length(data) > 0 &&
     all(vapply(data, function(z) inherits(z, cls), TRUE))) return(list(x=data, is_list=TRUE))
  stop("'data' must be ", what, " or a list of them (one per individual).")
}

# Individual ids for each row of each list element.
# Priority: 'id' column -> 'fallback' (e.g. ctmm identity) -> list names -> list index.
# A single (non-list) object without an 'id' column and without a fallback is individual "1".
proc_ids <- function(xl, id, is_list, fallback=NULL){
  nms <- names(xl)
  out <- vector("list", length(xl))
  for(k in seq_along(xl)){
    d <- xl[[k]]
    n <- nrow(d)
    if(!is.null(id) && id %in% names(d)){
      v <- as.character(d[[id]])
    } else {
      v <- if(!is.null(fallback)) fallback(d) else NULL
      if(is.null(v) || length(v) != 1 || is.na(v) || !nzchar(v)){
        v <- if(is_list){
          if(!is.null(nms) && !is.na(nms[k]) && nzchar(nms[k])) nms[k] else as.character(k)
        } else "1"
      }
      v <- rep(as.character(v), n)
    }
    out[[k]] <- v
  }
  out
}

# Build the 'times' data frame shared by proc_telem() and proc_Langevin().
# idv: character id per row; timestamp: stacked time stamps; tr: proc_resolve_time() result.
proc_times <- function(idv, timestamp, tr, max_gap){
  n <- length(idv)
  if(any(is.na(idv))) stop("Missing values in the individual ids.")
  if(any(duplicated(rle(idv)$values)))
    stop("Rows for each individual must be contiguous. Sort 'data' by id, then time ",
         "(or give each individual once in a list).")
  tnum <- as.numeric(timestamp)
  if(length(tnum) != n) stop("Time stamps are missing.")
  if(any(is.na(tnum))) stop("Missing values in the time stamps.")
  first <- !duplicated(idv)
  dts <- c(0, diff(tnum))
  dts[first] <- 0
  if(any(dts < 0)) stop("Observations are not in time order within individual.")
  dt <- dts * tr$scale
  new_seg <- first
  if(!is.null(max_gap)){
    if(!is.numeric(max_gap) || length(max_gap) != 1 || is.na(max_gap) || max_gap <= 0)
      stop("'max_gap' must be a single positive number (in 'time_unit').")
    new_seg <- new_seg | (dt > max_gap)
    dt[new_seg] <- 0
  }
  seg <- as.integer(cumsum(new_seg))
  times <- data.frame(obs=seq_len(n), id=idv, seg=seg, timestamp=timestamp, dt=dt,
                      stringsAsFactors=FALSE)
  attr(times, "time_unit") <- tr$unit
  attr(times, "ts_scale") <- tr$scale
  times
}

proc_check_return_type <- function(return_type){
  ok <- c("sparse", "dense", "data.frame", "matrix_df")
  if(!is.character(return_type) || length(return_type) != 1 || !(return_type %in% ok))
    stop("'return_type' must be one of ", paste0('"', ok, '"', collapse=", "), ".", call.=FALSE)
  invisible(NULL)
}

# Warn about nthreads when OpenMP is unavailable / limited.

proc_check_threads <- function(nthreads){
  if(nthreads > 1){
    oi <- omp_info()
    if(!oi$enabled){
      warning("'walk' was compiled without OpenMP support, so 'nthreads' is ignored and 1 thread is used. ",
              "See walk:::omp_info().", call.=FALSE)
    } else if(oi$thread_limit < nthreads){
      warning("OpenMP thread limit (OMP_THREAD_LIMIT) is ", oi$thread_limit, "; using at most that many threads.", call.=FALSE)
    }
  }
  invisible(NULL)
}

# Cell likelihoods + output assembly shared by proc_telem() and proc_Langevin().
# x, y: projected locations (units of cell_data); sxx, syy, sxy: error covariance;
# max_err: per-location neighbourhood radius (NA / <= 0 = no usable error information).
proc_lik <- function(x, y, sxx, syy, sxy, max_err, cell_data, times, return_type, trunc,
                     return_inadmissible, nthreads, neighborhood, clipRast){
  n <- length(x)
  valid <- !is.na(max_err) & max_err > 0
  me  <- ifelse(valid, max_err, 0)
  sxx <- ifelse(valid, sxx, 1)
  syy <- ifelse(valid, syy, 1)
  sxy <- ifelse(valid, sxy, 0)

  # Raster geometry and inadmissible-cell mask (NA in first layer)
  e <- as.vector(terra::ext(cell_data))
  geom <- c(e[["xmin"]], e[["ymax"]], terra::res(cell_data),
            terra::nrow(cell_data), terra::ncol(cell_data))
  mask <- is.na(terra::values(cell_data[[1]], mat=FALSE))

  # Optional clipping region: cells outside it are treated as inadmissible (like NA cells)
  win <- clip_window(clipRast, cell_data, x, y, me)
  if(!is.null(win)){
    mask <- clip_mask(mask, win)
    if(!win$data_box){
      n_cut <- clip_n_cut(win, cell_data, x, y, me)
      if(n_cut > 0) warning(n_cut, " location(s) have error neighbourhoods extending beyond 'clipRast'; ",
                            "likelihood mass outside it is dropped and the rest renormalised.", call.=FALSE)
    }
  }

  # Cell likelihoods (computed, masked, normalised and truncated in C++)
  if(neighborhood == "disc"){
    res <- bvn_lik_disc(x, y, sxx, syy, sxy, me, geom, mask, trunc, as.integer(nthreads))
  } else {
    vi <- which(valid)
    obs <- integer(0); cell <- numeric(0)
    if(length(vi)){
      t_buf <- terra::buffer(terra::vect(cbind(x[vi], y[vi]), crs=terra::crs(cell_data)), max_err[vi])
      extr <- terra::cells(cell_data, t_buf, touches = TRUE)
      obs  <- as.integer(vi[extr[, 'ID']])
      cell <- as.numeric(extr[, 'cell'])
      o <- order(obs); obs <- obs[o]; cell <- cell[o]
    }
    res <- bvn_lik_list(obs, cell, x, y, sxx, syy, sxy, me, geom, mask, trunc, as.integer(nthreads))
  }

  zero_ind <- length(res$inad_obs) > 0
  if(return_inadmissible){
    if(!zero_ind) return(NULL)
    return(cbind(obs=res$inad_obs, cell=res$inad_cell, lik=res$inad_lik))
  }
  if(zero_ind) warning("Some observations have likelihood values of 0 in all 'cell_data' cells! ",
                       "make_walk_data() excludes them (times$excluded == \"no_cells\").")

  ncells <- prod(dim(cell_data)[1:2])
  clip <- NULL
  if(!is.null(win)){
    res$cell <- clip_renumber(res$cell, win)
    ncells <- win$nrow * win$ncol
    clip <- clip_info(win)
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
  if(return_type == "data.frame"){
    out <- as.data.frame(out)
    out <- merge(out, times, by="obs")
    if(!is.null(clip)) attr(out, "clip") <- clip
    return(out)
  }
  out   # "matrix_df"
}

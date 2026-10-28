#' @title Predict Movement Path From Fitted CTMC Model
#' @param fit A fitted model object produced by \code{\link[walk]{fit_ctmc}}
#' @param walk_data Optional. The design data list (from \code{\link{make_walk_data}}) the
#' model was fitted to. Not needed, because \code{\link{fit_ctmc}} stores the observation
#' times in the fit; if supplied it is checked against the fit.
#' @param aux_timestamp Additional times for which location prediction is desired, of the same
#' kind as the data time stamps (POSIXct, or numeric for numeric `date`s). For POSIXct data,
#' `POSIXlt`/`Date`/character values are converted, numeric values are taken as seconds since
#' 1970-01-01 UTC, and all are shown in the time zone of the data. A single vector is used for
#' every segment (restricted to the segment's time span). For multi-individual data supply a
#' data frame with columns `id` and `timestamp` to give times per individual.
#' @param trunc_tol Level for truncating small location prediction probabilities. 
#' @param segments Optional integer vector of segment numbers (`walk_data$times$seg`) or
#' character (or factor) vector of individual ids to predict. Numeric values are always segment
#' numbers, so select numeric-looking ids with e.g. `segments = "1"`. Defaults to all.
#' @param debug Developer debugging. 
#' @details Each independent segment (individual, or burst created by `max_gap` in
#' \code{\link{proc_Langevin}}) is smoothed separately with the shared fitted rate matrix.
#' The smoothed state probabilities are returned as a sparse matrix (rows = prediction
#' times of all segments stacked, columns = cells) because a dense matrix is
#' `N x ns` and quickly exhausts memory.
#' @return A list (class `walk_pred`) with `local_state_prob` (sparse `dgCMatrix`), `times`, and
#' `cell_map`; use \code{\link{walk_rast}} to map rows to a raster. Column
#' `j` of `local_state_prob` is the cell with `cellx == j`; `cell_map$cell` gives its cell number in
#' the `cell_data` raster given to \code{\link{make_walk_data}} (also when
#' `clipRast` was used). Rows of `times` with `type == "p"` are `aux_timestamp` times; an
#' `aux_timestamp` equal to an observation time is not duplicated. `times$excluded` is the reason an
#' observation was excluded (`walk_data$times$excluded`, see \code{\link{make_walk_data}}; `NA` for
#' retained observations and `aux_timestamp` times). Predictions for excluded observations use the movement model only, like `aux_timestamp` times. As in
#' \code{\link{fit_ctmc}}, each segment starts at its first retained observation, so excluded
#' observations before it get the same distribution as that observation.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
predict_ctmc <- function(fit, walk_data=NULL, aux_timestamp=NULL, trunc_tol=1.0e-8, segments=NULL, debug=0){
  if(debug==1) browser()
  data_list <- fit$data_list
  if(is.null(walk_data)){
    if(is.null(fit$times)) stop("'fit' does not contain the observation times (it is not the output of ",
                                "fit_ctmc()); supply 'walk_data'.")
    walk_data <- list(times = fit$times, cell_data = fit$cell_data)
  }
  times_all <- walk_data$times
  if(is.null(times_all$seg)) times_all$seg <- 1L
  if(is.null(times_all$id)) times_all$id <- "1"
  time_unit <- attr(times_all, "time_unit")
  ts_scale <- proc_ts_scale(walk_data$times)
  times_all$row_obs <- seq_len(nrow(times_all))
  
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  par <- fit$par
  Xb_q_r <- data_list$X_q_r %*% par[par_map$beta_q_r]
  Xb_q_m <- if(length(par_map$beta_q_m)) data_list$X_q_m %*% par[par_map$beta_q_m] else matrix(0, length(data_list$from), 1)
  p <- if(!is.null(par_map$logit_p)) plogis(par[par_map$logit_p]) else 0
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
    delta <- get_lim_ud(fit)$ud; delta <- delta/sum(delta)
  } else delta <- data_list$delta
  Lt_all <- data_list$Lt
  if(is.null(Lt_all)) Lt_all <- Matrix::t(data_list$L)
  # observations impossible under the fitted model
  zz <- ctmc_n2ll_joint(par, data_list, by_segment = TRUE)
  zero_obs <- attr(zz, "zero_obs")
  if(length(zero_obs))
    stop(zero_obs_msg(zero_obs, times_all, "at the fitted parameters"), " No prediction is possible.", call. = FALSE)
  if(!all(is.finite(zz))) stop("The fitted model gives a non-finite likelihood; no prediction is possible.")
  
  # walk_data must be the one the model was fitted to
  if(nrow(times_all) != data_list$N || (!is.null(walk_data$q_r) && nrow(walk_data$q_r) != data_list$ns) ||
     (!is.null(walk_data$q_m) && nrow(walk_data$q_m) != length(data_list$from)) ||
     (!is.null(data_list$cell_map) && !is.null(walk_data$q_r) && !identical(as.numeric(walk_data$q_r$cell), as.numeric(data_list$cell_map$cell))))
    stop("'walk_data' does not match the data used for 'fit' (different observations or cells).")
  if(!is.null(fit$time_unit) && !identical(fit$time_unit, time_unit))
    stop("'walk_data' time_unit (", time_unit, ") differs from the fitted model (", fit$time_unit, ").")
  
  segs <- unique(times_all$seg)
  if(!is.null(segments)){
    if(is.character(segments) || is.factor(segments)){
      segments <- as.character(segments)
      miss <- setdiff(segments, times_all$id)
      if(length(miss)) warning("'segments' ids not found in 'walk_data$times$id': ", paste(miss, collapse=", "))
      segs <- unique(times_all$seg[times_all$id %in% segments])
    } else {
      miss <- setdiff(segments, segs)
      if(length(miss)) warning("'segments' numbers not found in 'walk_data$times$seg': ", paste(miss, collapse=", "),
                               ". Use a character vector to select individuals by id.")
      segs <- intersect(segs, segments)
    }
    if(!length(segs)) stop("No segments selected.")
  }
  # aux_timestamp: same class / time zone as the data time stamps
  ts0 <- times_all$timestamp
  as_ts <- function(a){
    if(inherits(ts0, "POSIXt")){
      if(inherits(a, "POSIXt")) a <- as.POSIXct(a)
      else if(inherits(a, "Date")) a <- as.POSIXct(a)
      else if(is.numeric(a)) a <- as.POSIXct(a, origin = "1970-01-01")
      else if(is.character(a)) a <- as.POSIXct(a, tz = if(is.null(attr(ts0, "tzone"))) "" else attr(ts0, "tzone"))
      else stop("'aux_timestamp' must be POSIXct (like the data time stamps).")
      attr(a, "tzone") <- attr(ts0, "tzone")
    } else if(is.numeric(ts0)){
      if(!is.numeric(a)) stop("The data time stamps are numeric, so 'aux_timestamp' must be numeric too.")
      a <- as.numeric(a)
    }
    a
  }
  if(!is.null(aux_timestamp)){
    if(is.data.frame(aux_timestamp)){
      if(!all(c("id","timestamp") %in% names(aux_timestamp))) stop("An 'aux_timestamp' data frame needs columns 'id' and 'timestamp'.")
      aux_timestamp <- data.frame(id = as.character(aux_timestamp$id), stringsAsFactors = FALSE,
                                  timestamp = as_ts(aux_timestamp$timestamp))
    } else aux_timestamp <- as_ts(aux_timestamp)
  }
  if(debug==2) browser()
  
  res_P <- vector("list", length(segs)); res_T <- vector("list", length(segs))
  aux_used <- character(0)
  for(k in seq_along(segs)){
    tk <- times_all[times_all$seg == segs[k], c("obs","id","seg","timestamp","dt","row_obs")]
    tk$type <- "o"
    # reason for exclusion (NA = retained); rows of L without information are excluded
    ex <- if(is.null(times_all$excluded)) rep(NA_character_, nrow(tk)) else as.character(times_all$excluded[tk$row_obs])
    ex[is.na(ex) & Matrix::colSums(Lt_all[, tk$row_obs, drop=FALSE]) == 0] <- "no_cells"
    tk$excluded <- ex
    if(!is.null(aux_timestamp)){
      if(is.data.frame(aux_timestamp)){
        ak <- aux_timestamp$timestamp[aux_timestamp$id == tk$id[1]]
      } else ak <- aux_timestamp
      ak <- ak[ak >= min(tk$timestamp) & ak <= max(tk$timestamp)]
      aux_used <- c(aux_used, paste(tk$id[1], as.numeric(ak)))
      if(length(ak)){
        tk <- rbind(tk, data.frame(obs=NA, id=tk$id[1], seg=tk$seg[1], timestamp=ak, dt=NA, row_obs=NA, type="p", excluded=NA_character_))
      }
    }
    tk <- tk[order(tk$timestamp, tk$type),]
    # an aux time equal to an observation time is that observation (no extra zero-dt row)
    tk <- tk[!(tk$type == "p" & (tk$timestamp %in% tk$timestamp[tk$type == "o"] | duplicated(tk$timestamp))),]
    tk$dt <- c(0, diff(as.numeric(tk$timestamp))) * ts_scale
    # as in fit_ctmc(), the segment starts at its first retained observation
    fr <- which(tk$type == "o" & is.na(tk$excluded))
    tk$dt[seq_len(if(length(fr)) fr[1] else nrow(tk))] <- 0
    Lk <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(data_list$ns, nrow(tk)))
    io <- which(tk$type=="o")
    Lk[, io] <- Lt_all[, tk$row_obs[io], drop=FALSE]
    G <- ctmc_predict_arma(
      L=as(Lk, "dgCMatrix"), dt=tk$dt, ns=data_list$ns, from_to=from_to,
      Xb_q_r=Xb_q_r, Xb_q_m=Xb_q_m, p=p, delta=matrix(delta, nrow=1),
      hij=data_list$hij, eq_prec=data_list$eq_prec, trunc_tol=trunc_tol,
      link_r=which(data_list$link_r==c("soft_plus","log","logit")),
      a_r=data_list$a_r, l_r=data_list$l_r, u_r=data_list$u_r,
      link_m=which(data_list$link_m==c("soft_plus","log")),
      a_m=data_list$a_m, form=which(data_list$form==c("mult","sde")),
      k=data_list$k, norm=data_list$norm, clip=data_list$clip)
    # truncate small values and renormalise so that every row is a probability distribution
    G[G < trunc_tol] <- 0
    rs <- rowSums(G); G <- G / ifelse(rs > 0, rs, 1)
    res_P[[k]] <- as(Matrix::Matrix(G, sparse=TRUE), "dgCMatrix")
    rm(G)
    res_T[[k]] <- tk
  }
  # aux times not inside any selected segment (before the first / after the last fix of an
  # individual, or inside a 'max_gap' gap) are not predicted: say so
  if(!is.null(aux_timestamp)){
    ids_sel <- unique(times_all$id[times_all$seg %in% segs])
    if(is.data.frame(aux_timestamp)){
      ak_all <- paste(aux_timestamp$id, as.numeric(aux_timestamp$timestamp))[aux_timestamp$id %in% ids_sel]
    } else ak_all <- as.vector(outer(ids_sel, as.numeric(aux_timestamp), paste))
    n_drop <- sum(!(unique(ak_all) %in% aux_used))
    if(n_drop) warning(n_drop, " 'aux_timestamp' time(s) (per id) lie outside the time span of every selected ",
                       "segment (before the first / after the last observation, or in a 'max_gap' gap) and were not predicted.",
                       call.=FALSE)
  }
  times <- do.call(rbind, res_T)
  times$row <- seq_len(nrow(times))
  attr(times, "time_unit") <- time_unit
  attr(times, "ts_scale") <- ts_scale
  if(debug==3) browser()
  out <- list(local_state_prob=do.call(rbind, res_P), times=times, cell_map=data_list$cell_map)
  attr(out, "cell_data") <- if(!is.null(fit$cell_data)) fit$cell_data else walk_data$cell_data
  class(out) <- c("walk_pred", "list")
  return(out)
}

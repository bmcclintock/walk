#' @title Predict Movement Path From Fitted CTMC Model
#' @param fit A fitted model object produced by \code{\link[walk]{fit_ctmc}}
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param aux_timestamp Additional times for which location prediction is desired, of the same
#' kind as the data time stamps (POSIXct, or numeric for numeric `date`s). For
#' multi-individual data supply a data frame with columns `id` and `timestamp`.
#' @param trunc_tol Level for truncating small location prediction probabilities. 
#' @param segments Optional integer vector of segment numbers (`walk_data$times$seg`) or
#' character vector of individual ids to predict. Defaults to all.
#' @param debug Developer debugging. 
#' @param ... Aditional arguments passed to internal function \code{ctmc_predict_arma}. 
#' @details Each independent segment (individual, or burst created by `max_gap` in
#' \code{\link{proc_Langevin}}) is smoothed separately with the shared fitted rate matrix.
#' The smoothed state probabilities are returned as a sparse matrix (rows = prediction
#' times of all segments stacked, columns = cells) because a dense matrix is
#' `N x ns` and quickly exhausts memory.
#' @return A list with `local_state_prob` (sparse `dgCMatrix`) and `times`.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
predict_ctmc <- function(fit, walk_data, aux_timestamp=NULL, trunc_tol=1.0e-8, segments=NULL, debug=0, ...){
  if(debug==1) browser()
  data_list <- fit$data_list
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
  
  segs <- unique(times_all$seg)
  if(!is.null(segments)){
    if(is.character(segments)) segs <- unique(times_all$seg[times_all$id %in% segments])
    else segs <- intersect(segs, segments)
  }
  if(debug==2) browser()
  
  res_P <- vector("list", length(segs)); res_T <- vector("list", length(segs))
  for(k in seq_along(segs)){
    tk <- times_all[times_all$seg == segs[k], c("obs","id","seg","timestamp","dt","row_obs")]
    tk$type <- "o"
    if(!is.null(aux_timestamp)){
      if(is.data.frame(aux_timestamp)){
        ak <- aux_timestamp$timestamp[aux_timestamp$id == tk$id[1]]
      } else ak <- aux_timestamp
      ak <- ak[ak >= min(tk$timestamp) & ak <= max(tk$timestamp)]
      if(length(ak)){
        tk <- rbind(tk, data.frame(obs=NA, id=tk$id[1], seg=tk$seg[1], timestamp=ak, dt=NA, row_obs=NA, type="p"))
      }
    }
    tk <- tk[order(tk$timestamp, tk$type),]
    tk <- tk[!duplicated(tk[,c("timestamp","type")]),]
    tk$dt <- c(0, diff(as.numeric(tk$timestamp))) * ts_scale
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
    G[G < trunc_tol] <- 0
    res_P[[k]] <- as(Matrix::Matrix(G, sparse=TRUE), "dgCMatrix")
    rm(G)
    res_T[[k]] <- tk
  }
  times <- do.call(rbind, res_T)
  times$row <- seq_len(nrow(times))
  attr(times, "time_unit") <- time_unit
  attr(times, "ts_scale") <- ts_scale
  if(debug==3) browser()
  return(list(local_state_prob=do.call(rbind, res_P), times=times))
}

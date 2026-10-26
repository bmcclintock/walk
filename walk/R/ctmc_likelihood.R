#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate plogis
#' @importFrom Matrix diag
#' @author Devin S. Johnson
#' @export
ctmc_n2ll <- function(par, data_list, debug=0, ...){
  if(debug==1) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  
  if(!is.null(par_map$logit_p)){
    p = plogis(logit_p)
  } else{
    p = 0
  }
  
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
    delta <- get_lim_ud(list(par = par, data_list = data_list))
    delta <- delta$ud
    delta <- delta/sum(delta)
  } else{
    delta <- data_list$delta
  }
  
  if(debug==2) browser()
  
  ctmc_n2ll_arma(
    L = t(data_list$L), #arma code needs the obs on the columns
    dt = data_list$dt, 
    ns = data_list$ns, 
    from_to = from_to, 
    Xb_q_r = Xb_q_r, 
    Xb_q_m = Xb_q_m,
    p = p,
    delta = matrix(delta, nrow=1),
    hij = data_list$h,
    eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
    a_r = data_list$a_r,
    l_r = data_list$l_r,
    u_r = data_list$u_r,
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    a_m = data_list$a_m,
    form = which(data_list$form==c("mult", "sde")),
    k = data_list$k,
    norm = data_list$norm,
    clip = data_list$clip
  )
  
}




#' @title Joint (multi-individual) movement MMPP log-likelihood
#' @description Evaluates the CTMC likelihood for data containing several
#' independent chains (individuals, or segments of individuals split at long gaps)
#' that share one set of movement parameters. The total is the sum of the
#' per-segment log-likelihoods, each segment starting from `delta`.
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood. Must contain
#' `seg_start` (0-based first observation of each segment). If absent the data are
#' treated as one segment.
#' @param by_segment Logical. Return the vector of per-segment -2 log-likelihoods
#' (with attribute `zero_obs`) instead of their sum.
#' @param debug For developers only, leave in the default setting.
#' @param ... Ignored.
#' @details The forward recursion only touches cells that can carry probability
#' (the support of the observation likelihood, expanded by the uniformisation
#' steps), which is exact, and segments are evaluated in parallel with OpenMP using
#' `data_list$nthreads` threads. If `data_list$skip_zero` is `TRUE`, observations
#' whose error support does not overlap the predicted distribution are treated as
#' missing instead of returning `Inf`; their indices are returned in
#' `attr(, "zero_obs")`.
#' @author Devin S. Johnson, Brett T. McClintock
#' @export
ctmc_n2ll_joint <- function(par, data_list, by_segment=FALSE, debug=0, ...){
  if(debug==1) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  logit_p <- par[par_map$logit_p]
  
  Xb_q_r <- as.vector(data_list$X_q_r %*% beta_q_r)
  Xb_q_m <- if(length(beta_q_m)) as.vector(data_list$X_q_m %*% beta_q_m) else rep(0, length(data_list$from))
  p <- if(!is.null(par_map$logit_p)) plogis(logit_p) else 0
  
  if(is.character(data_list$delta) && data_list$delta=="stationary"){
    delta <- get_lim_ud(list(par = par, data_list = data_list))$ud
    delta <- delta/sum(delta)
  } else{
    delta <- data_list$delta
  }
  seg_start <- data_list$seg_start
  if(is.null(seg_start)) seg_start <- 0L
  Lt <- data_list$Lt
  if(is.null(Lt)) Lt <- Matrix::t(data_list$L)
  nthreads <- data_list$nthreads; if(is.null(nthreads)) nthreads <- 1L
  prune <- data_list$prune; if(is.null(prune)) prune <- 0
  skip_zero <- isTRUE(data_list$skip_zero)
  if(debug==2) browser()
  r <- ctmc_n2ll_joint_arma(
    L = Lt, dt = data_list$dt, seg_start = as.integer(seg_start),
    ns = data_list$ns, from_to = from_to,
    Xb_q_r = Xb_q_r, Xb_q_m = Xb_q_m, p = p,
    delta = matrix(delta, nrow=1), hij = data_list$hij,
    eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
    a_r = data_list$a_r, l_r = data_list$l_r, u_r = data_list$u_r,
    link_m = which(data_list$link_m==c("soft_plus", "log")),
    a_m = data_list$a_m,
    form = which(data_list$form==c("mult", "sde")),
    k = data_list$k, norm = data_list$norm, clip = data_list$clip,
    nthreads = as.integer(nthreads), prune = prune, skip_zero = skip_zero
  )
  out <- if(by_segment) r$n2ll else sum(r$n2ll)
  if(length(r$zero_obs)) attr(out, "zero_obs") <- r$zero_obs
  out
}

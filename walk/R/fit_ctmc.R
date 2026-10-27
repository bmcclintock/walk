#' @title Fit CTMC movement model to telemetry data
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param model_parameters A named list giving specifications for the model. See \code{\link[walk]{ctmc_control}}. 
#' @param pen_fun An optional penalty function. Should be on the scale of a log-prior distribution.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q_r=c(), beta_q_r=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param eq_prec Error rate of matrix exponential calculation. Defaults to \code{1.0e-8}. This is 
#' a generous value. If the model is running slow, you can try reducing it to, say, \code{1.0e-4}.
#' @param nthreads Number of OpenMP threads used to evaluate independent segments
#' (individuals / bursts) of the joint likelihood in parallel. Defaults to 1.
#' @param prune Optional threshold for dropping filtered state probabilities smaller than
#' `prune * max(phi)` during the forward recursion. `0` (default) is exact.
#' @param skip_zero Logical. If `TRUE`, observations whose error support has zero predictive
#' probability are treated as missing rather than making the likelihood infinite. The indices
#' of such observations at the solution are returned in `zero_obs`.
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @details Data from several individuals (see \code{\link{proc_Langevin}}, argument `id`, and
#' \code{\link{bind_proc_data}}) are fitted as a single joint model: all individuals share
#' the movement parameters, each individual (or each segment created with `max_gap`) is an
#' independent chain started from `delta`, and the objective is the sum of the per-segment
#' -2 log-likelihoods. A single-individual data set is the special case with one segment
#' and gives exactly the same value as before. Rates are per unit of `walk_data$times$dt`,
#' i.e. per `time_unit` of \code{\link{proc_Langevin}} / \code{\link{proc_telem}}; the unit is
#' returned as `time_unit`.
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
fit_ctmc <- function(walk_data, 
                     model_parameters = ctmc_control(), 
                     pen_fun = NULL, hessian=TRUE, start=NULL, method="nlminb", 
                     fit=TRUE, eq_prec = 1.0e-8, nthreads = 1, prune = 0,
                     skip_zero = FALSE, debug=0, ...){
  
  if(debug==1) browser()
  # cell_idx_df <- select(walk_data$q_r, cell, cellx) %>% distinct()
  # data <- data %>% left_join(cell_idx_df, by="cell")
  
  if(is.null(model_parameters$q_r$form) | is.null(model_parameters$q_m$form)) stop("Model formulas must be specified for both residency and movement.")
  X_q_r <- dm_q_r(model_parameters$q_r$form, walk_data)
  X_q_m <- dm_q_m(model_parameters$q_m$form, walk_data)
  
  par_map <- list(beta_q_r = 1:ncol(X_q_r))
  if(ncol(X_q_m)!=0) par_map$beta_q_m <- c(1:ncol(X_q_m)) + ncol(X_q_r)
  
  if(is.null(model_parameters$p)) model_parameters$p <- FALSE
  if(model_parameters$p){
    par_map$logit_p <- ncol(X_q_m) + ncol(X_q_r) + 1
  }
  
  delta <- model_parameters$delta
  if(is.numeric(delta)){
    if(length(delta) != nrow(walk_data$q_r)) stop("The length of 'delta' vector is not equal to the number of cells.")
    delta <- delta/sum(delta)
  }
  if(delta=="uniform"){
    delta <- rep(1,nrow(walk_data$q_r))
    delta <- delta/sum(delta)
  }
  
  if(is.null(model_parameters$q_r$link)){
    link_r <- "soft_plus"
  } else{
    link_r <- model_parameters$q_r$link
  }
  if(is.null(model_parameters$q_m$link)){
    link_m <- "soft_plus"
  } else{
    link_m <- model_parameters$q_m$link
  }
  if(!link_r%in%c("soft_plus","log","logit")) stop("The 'link_r' objects in must be either 'soft_plus', 'log', or 'logit'")
  if(!link_m%in%c("soft_plus","log")) stop("The 'link_m' objects in must be either 'soft_plus' or 'log'.")
  
  
  form <- model_parameters$form
  if(!is.null(form)){
    if(!form%in%c("mult","sde")) stop("The 'form' object in must be either 'mult',or 'sde'.")
  } else{
    form<- "mult"
  }
  
  
  if(is.null(model_parameters$q_r$a)){
    a_r <- 1
  } else{
    a_r <- model_parameters$q_r$a
  }
  
  if(is.null(model_parameters$q_r$U)){
    u_r <- 0
  } else{
    u_r <- model_parameters$q_r$U
  }
  if(is.null(model_parameters$q_r$L)){
    l_r <- 0
  } else{
    l_r <- model_parameters$q_r$L
  }
  if(link_r=="logit"){
    if(l_r < 0) stop("L for the logit link must be >= 0")
    if(u_r <= 0) stop("U for the logit link must be > 0")
    if(u_r <= l_r) stop("U must be > L for the logit link")
  }
  if(is.null(model_parameters$q_m$a)){
    a_m <- 1
  } else{
    a_m <- model_parameters$q_m$a
  }
  
  if(a_r<1 | a_m<1) stop("The 'a' parameter for the 'soft_plus' and 'sde' link functions must be >1.")
  
  if(form %in% c("mult","add")){
    k <- 0
  } else {
    k <- if(attr(walk_data$q_m, "directions")=="rook"){
      k <- 2
    } else if(attr(walk_data$q_m, "directions")=="queen"){
      k <- 4
    } else {
      stop("'sde' model currently only functions with 'SpatRaster' class habitat data.")
    }
  }
  
  norm <- model_parameters$norm
  if(is.null(norm) | !is.logical(norm)) norm <- TRUE
  if(model_parameters$form=="sde") norm <- FALSE
  
  clip <- model_parameters$clip
  if(is.null(clip) || clip<0) clip <- 0
  
  times <- walk_data$times
  seg <- times$seg
  if(is.null(seg)) seg <- rep(1L, nrow(times))
  seg_start <- which(c(TRUE, seg[-1] != seg[-length(seg)])) - 1L
  if(any(duplicated(seg[seg_start + 1L]))) stop("Segments in 'walk_data$times' are not contiguous.")
  dt <- times$dt
  dt[seg_start + 1L] <- 0
  if(any(dt < 0)) stop("Negative time differences within a segment. Re-run proc_telem() or proc_Langevin() with the 'id' argument.")
  time_unit <- attr(times, "time_unit")
  
  data_list <- list(
    N = nrow(walk_data$L),
    seg_start = as.integer(seg_start),
    Lt = as(Matrix::t(walk_data$L), "dgCMatrix"),
    nthreads = as.integer(nthreads),
    prune = prune,
    skip_zero = skip_zero,
    ns = nrow(walk_data$q_r),
    dt = dt,
    L = walk_data$L,
    delta = delta,
    hij=walk_data$q_m$h,
    ### Q
    from = as.integer(walk_data$q_m$from_cellx-1),
    to = as.integer(walk_data$q_m$cellx-1),
    X_q_r = X_q_r,
    X_q_m = X_q_m,
    par_map = par_map,
    eq_prec = eq_prec,
    link_r=link_r,
    a_r = a_r,
    l_r = l_r,
    u_r = u_r,
    link_m=link_m,
    a_m = a_m,
    form=form,
    k = k,
    norm = norm,
    clip = clip,
    cell_map = walk_data$q_r[,c("cell","cellx")]
  )
  
  
  if(is.null(start$beta_q_r)) start$beta_q_r <- rep(0, ncol(X_q_r))
  if(is.null(start$beta_q_m)) start$beta_q_m <- rep(0, ncol(X_q_m))
  if(is.null(start$logit_p) & model_parameters$p) start$logit_p <- qlogis(0.05)
  
  par_start <- c(start$beta_q_r, start$beta_q_m, start$logit_p)
  
  # ctmc_n2ll(par, data_list)
  
  if(is.null(pen_fun)){
    obj_fun <- function(par, data_list, ...){as.numeric(ctmc_n2ll_joint(par, data_list, debug=0))}
  } else{
    obj_fun <- function(par, data_list, ...){as.numeric(ctmc_n2ll_joint(par, data_list, debug=0)) - 2*pen_fun(par)}
  }
  
  if(debug==2) browser()
  
  # opt <- optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, control = list(trace=10, rel.tol=1.0e-3))
  # ctmc_n2ll(par_start, data_list, debug=1)
  
  if(fit){
    message('Optimizing likelihood...')  
    if(debug==2) browser()
    opt <- optimx::optimr(par=par_start, fn=obj_fun, method=method, data_list=data_list, ...)
    
    if(opt$convergence!=0){
      warning("There might be a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      # hessian <- FALSE
      # V <- NULL
    }
    if(hessian){
      message('Calculating Hessian and variance-covariance matrices...')  
      H <- numDeriv::hessian(obj_fun, opt$par, data_list=data_list)
      V <- tryCatch(2*solve(H), error=function(e){
        warning("Hessian is singular; standard errors not available.")
        NULL})
    } else{
      V <- NULL
      H <- NULL 
    }
  } else{
    # do not evaluate the (possibly very expensive) likelihood when fit = FALSE
    hessian <- FALSE
    V <- NULL
    opt <- list(par=par_start, value=NA_real_, convergence=NA_integer_)
  }
  
  if(debug==3) browser()
  
  par <- as.vector(opt$par)
  beta <- get_betas(par, V, data_list)
  
  if(!hessian){
    V <- NULL
    H <- NULL
  }
  
  zero_obs <- NULL
  if(fit){
    zz <- ctmc_n2ll_joint(par, data_list, by_segment=TRUE)
    zero_obs <- attr(zz, "zero_obs")
    seg_n2ll <- as.vector(zz)
  } else seg_n2ll <- NULL
  
  out <- list(
    par = par,
    seg_n2ll = seg_n2ll,
    zero_obs = zero_obs,
    n_seg = length(seg_start),
    time_unit = time_unit,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    hessian = 0.5*H,
    results = beta, #list(beta = beta, real = reals),
    opt = opt,
    start = start,
    data_list = data_list,
    pen_fun = pen_fun,
    model_parameters = model_parameters
  )
  
  if(debug==4) browser()
  
  return(out)
  
}
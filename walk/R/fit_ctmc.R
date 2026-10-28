#' @title Fit CTMC movement model to telemetry data
#' @param walk_data A design data list produced by the function \code{\link{make_walk_data}}.
#' @param model_parameters A named list giving specifications for the model. See \code{\link[walk]{ctmc_control}}. 
#' @param pen_fun An optional penalty function. Should be on the scale of a log-prior distribution.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values: a named list with elements `beta_q_r`, `beta_q_m` and,
#' if `p = TRUE`, `logit_p` (any may be omitted; omitted ones start at 0, or `qlogis(0.05)` for
#' `logit_p`). Each must have one value per model-matrix column, in the order of the columns
#' (e.g. `(Intercept)` first). A fitted model's `fit$start` has this form.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param eq_prec Error rate of matrix exponential calculation. Defaults to \code{1.0e-8}. This is 
#' a generous value. If the model is running slow, you can try reducing it to, say, \code{1.0e-4}.
#' @param nthreads Number of OpenMP threads used to evaluate independent segments
#' (individuals / bursts) of the joint likelihood in parallel. Defaults to 1.
#' @param prune Optional threshold for dropping filtered state probabilities smaller than
#' `prune * max(phi)` during the forward recursion. `0` (default) is exact.
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @details Data from several individuals (see \code{\link{proc_Langevin}}, argument `id`, and
#' \code{\link{bind_proc_data}}) are fitted as a single joint model: all individuals share
#' the movement parameters, each individual (or each segment created with `max_gap`) is an
#' independent chain started from `delta`, and the objective is the sum of the per-segment
#' -2 log-likelihoods. A single-individual data set is the special case with one segment.
#' Rates are per unit of `walk_data$times$dt`,
#' i.e. per `time_unit` of \code{\link{proc_Langevin}} / \code{\link{proc_telem}}; the unit is
#' returned as `time_unit`. `seg_info` gives the segment number, individual `id`, first row
#' and number of observations of each segment (and its -2 log-likelihood `n2ll` when
#' `fit = TRUE`); `seg_n2ll` is named by segment number. `n_excluded` counts the excluded
#' observations (`walk_data$times$excluded`, set by \code{\link{make_walk_data}} and
#' \code{\link{exclude_obs}}). An excluded observation is exactly equivalent to a deleted one: the
#' next time step is measured from the previous retained observation, and a segment starts (from
#' `delta`) at its first retained observation.
#'
#' An observation that cannot be reached from the previous ones under the model (e.g. a gross
#' outlier) gives likelihood zero. If this happens at the starting values, `fit_ctmc()` stops
#' with an error that lists the first such observation of each affected segment; exclude them with
#' `make_walk_data(vmax = )` or \code{\link{exclude_obs}} (or, if they are plausible, use starting
#' values with a higher rate). `zero_obs` gives such observations at the returned parameters.
#' @param ... Additional named arguments passed to the optimization function
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}: `control`,
#' `lower`, `upper`, `gr` or `hess`. Any other argument is an error.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois qlogis
#' @export
fit_ctmc <- function(walk_data, 
                     model_parameters = ctmc_control(), 
                     pen_fun = NULL, hessian=TRUE, start=NULL, method="nlminb", 
                     fit=TRUE, eq_prec = 1.0e-8, nthreads = 1, prune = 0,
                     debug=0, ...){
  
  if(debug==1) browser()
  # '...' goes to optimx::optimr(), which does not reject unknown arguments (they reach the
  # objective and can make the optimiser fail silently), so only its own arguments are allowed
  ok_dots <- c("gr", "hess", "lower", "upper", "control")
  nm <- names(list(...)); if(is.null(nm)) nm <- rep("", ...length())
  if(any(!(nm %in% ok_dots)))
    stop("Unused argument(s) in fit_ctmc(): ",
         paste(ifelse(nzchar(nm), nm, "<unnamed>")[!(nm %in% ok_dots)], collapse=", "),
         ". Arguments in '...' are passed to optimx::optimr() and must be one of: ",
         paste(ok_dots, collapse=", "), ".", call. = FALSE)
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
  if(is.null(delta)) delta <- "uniform"
  if(is.numeric(delta)){
    if(length(delta) != nrow(walk_data$q_r)) stop("The length of 'delta' vector is not equal to the number of cells (nrow(walk_data$q_r)).")
    if(any(!is.finite(delta)) || any(delta < 0) || sum(delta) <= 0) stop("A numeric 'delta' must be non-negative and finite with a positive sum.")
    delta <- delta/sum(delta)
  } else if(identical(delta, "uniform")){
    delta <- rep(1,nrow(walk_data$q_r))
    delta <- delta/sum(delta)
  } else if(!identical(delta, "stationary")){
    stop("'delta' must be \"uniform\", \"stationary\", or a numeric vector with one value per cell.")
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
  # excluded observations (empty rows of L); a segment starts at its first retained observation,
  # so excluded observations before it are the same as deleted ones
  excl <- Matrix::rowSums(walk_data$L) == 0
  dt[seg_lead(excl, seg)] <- 0
  time_unit <- attr(times, "time_unit")
  
  data_list <- list(
    N = nrow(walk_data$L),
    seg_start = as.integer(seg_start),
    Lt = as(Matrix::t(walk_data$L), "dgCMatrix"),
    nthreads = as.integer(nthreads),
    prune = prune,
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
  
  
  if(!is.null(start)){
    if(!is.list(start) || is.null(names(start)) || any(!nzchar(names(start))))
      stop("'start' must be a named list with elements 'beta_q_r', 'beta_q_m' (and 'logit_p' if p = TRUE).")
    bad <- setdiff(names(start), c("beta_q_r","beta_q_m","logit_p"))
    if(length(bad)) stop("Unknown element(s) in 'start': ", paste(bad, collapse=", "),
                         ". Use 'beta_q_r', 'beta_q_m' and 'logit_p'.")
    nb <- c(beta_q_r=ncol(X_q_r), beta_q_m=ncol(X_q_m), logit_p=as.integer(isTRUE(model_parameters$p)))
    for(nm in intersect(names(start), names(nb))){
      if(length(start[[nm]]) != nb[[nm]]) stop("'start$", nm, "' has length ", length(start[[nm]]), " but the model has ",
                                                 nb[[nm]], " parameter(s): ",
                                                 paste(switch(nm, beta_q_r=colnames(X_q_r), beta_q_m=colnames(X_q_m), logit_p="logit_p"), collapse=", "), ".")
    }
  }
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
    z0 <- ctmc_n2ll_joint(par_start, data_list, by_segment=TRUE)
    if(!all(is.finite(z0))) stop(zero_obs_msg(attr(z0, "zero_obs"), times, "at the starting values"), call. = FALSE)
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
  seg_info <- data.frame(seg = seg[seg_start + 1L],
                         id = if(is.null(times$id)) "1" else as.character(times$id[seg_start + 1L]),
                         start = seg_start + 1L, n = diff(c(seg_start, length(seg))),
                         stringsAsFactors = FALSE)
  # excluded observations (see make_walk_data(), exclude_obs())
  seg_info$n_excluded <- as.vector(tapply(excl, factor(seg, levels = seg_info$seg), sum))
  if(!is.null(seg_n2ll)){
    seg_info$n2ll <- seg_n2ll
    names(seg_n2ll) <- seg_info$seg
  }
  
  out <- list(
    par = par,
    seg_n2ll = seg_n2ll,
    zero_obs = zero_obs,
    n_seg = length(seg_start),
    seg_info = seg_info,
    time_unit = time_unit,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    hessian = 0.5*H,
    results = beta, #list(beta = beta, real = reals),
    opt = opt,
    start = start,
    data_list = data_list,
    times = walk_data$times,          # so predict_ctmc() does not need 'walk_data'
    cell_data = walk_data$cell_data,  # raster template (wrapped) for walk_rast()
    pen_fun = pen_fun,
    model_parameters = model_parameters
  )
  
  if(debug==4) browser()
  
  return(out)
  
}
# TRUE for the rows of each segment up to and including its first retained (non-excluded) row
#' @noRd
seg_lead <- function(excl, seg){
  first <- c(TRUE, seg[-1] != seg[-length(seg)])
  n_kept <- stats::ave(as.integer(!excl), cumsum(first), FUN = cumsum)
  n_kept <= 1 & !(n_kept == 1 & excl)
}

# error message for observations with zero likelihood ('zero_obs' rows of 'times')
#' @noRd
zero_obs_msg <- function(zero_obs, times, where){
  if(!length(zero_obs)) return(paste0("The likelihood is not finite ", where, "."))
  id <- if(is.null(times$id)) "1" else as.character(times$id[zero_obs])
  lab <- paste0(zero_obs, " (", id, ")")
  paste0("The likelihood is zero ", where, ": observation(s) ",
         paste(utils::head(lab, 10), collapse=", "), if(length(lab) > 10) ", ...",
         " (rows of walk_data$times; id) cannot be reached from the preceding observations of their ",
         "segment. Only the first such observation of each segment is listed. These are usually outliers: ",
         "exclude them with make_walk_data(..., vmax = <maximum speed>) or ",
         "exclude_obs(walk_data, c(", paste(utils::head(zero_obs, 10), collapse=", "),
         if(length(zero_obs) > 10) ", ..." , ")) and refit. If they are plausible, try starting values with a ",
         "higher movement rate ('start$beta_q_r').")
}

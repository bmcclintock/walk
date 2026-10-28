#' @title Create sparse representations of design matrices for 
#' movement and resight models
#' @param formula An \code{R} formula object for the paramter vector.
#' @param ddl Design data list
#' @details This function is not designed for end-users but is exported for posterity.
#' @name dm_matrix
NULL

#' @rdname dm_matrix
#' @importFrom stats model.matrix sd
#' @export
dm_q_m <- function(formula, ddl){
  X <- dm_check_na(formula, ddl$q_m, "q_m")
  # if(is.null(par_list$offset)){
  #   par_list$offset <- ~0 + log(1/num_neigh)
  #   offset = model.matrix(par_list$offset, ddl$q)
  #   offset = as.vector(offset)
  # } else if(par_list$offset==0){
  #   offset=rep(0,nrow(X))
  # }else{
  #   offset = model.matrix(par_list$offset, ddl$q)
  #   if(ncol(offset)>1) stop("q offset formula must result in model.matrix of only 1 column.")
  #   offset = as.vector(offset)
  #   offset = ifelse(is.na(offset), 0, offset)
  # }
  keep_col <- !(apply(X, 2, sd)==0)
  X <- X[,keep_col,drop=FALSE]
  #if(ncol(X)==0) X <- NA
  # uX <- unique(X)
  # uX <- uX[rowSums(uX)!=0,,drop=FALSE]
  # dX <- data.frame(cbind(from_cellx=ddl$q$from_cellx, to_cellx=ddl$q$to_cellx, X))
  # duX <- data.frame(cbind(idx_q=1:nrow(uX), uX))
  # mX <- merge(dX, duX)
  # mX <- with(mX, mX[order(from_cellx, to_cellx),])
  # lookup <- mX[,c('from_cellx','to_cellx','idx_q')]
  # return(list(X_q = uX, idx_q=lookup, off_q=offset))
  return(X)
}

#' @rdname dm_matrix
#' @export
dm_q_r <- function(formula, ddl){
  X <- dm_check_na(formula, ddl$q_r, "q_r")
  keep_col <- !((colMeans(X)!=1) & (apply(X, 2, sd)==0))
  X <- X[,keep_col,drop=FALSE]
  
  #if(ncol(X)==0) X <- NA
  # uX <- unique(X)
  # uX <- uX[rowSums(uX)!=0,,drop=FALSE]
  # dX <- data.frame(cbind(from_cellx=ddl$q$from_cellx, to_cellx=ddl$q$to_cellx, X))
  # duX <- data.frame(cbind(idx_q=1:nrow(uX), uX))
  # mX <- merge(dX, duX)
  # mX <- with(mX, mX[order(from_cellx, to_cellx),])
  # lookup <- mX[,c('from_cellx','to_cellx','idx_q')]
  # return(list(X_q = uX, idx_q=lookup, off_q=offset))
  return(X)
}

# model.matrix() silently drops rows with NA covariates, which would misalign the design matrix
# with the cells (q_r) or cell pairs (q_m). Keep all rows and stop with an informative error.
#' @noRd
dm_check_na <- function(formula, data, what){
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  X <- model.matrix(formula, mf)
  bad <- !stats::complete.cases(X)
  if(any(bad)){
    v <- all.vars(formula)
    v <- v[vapply(v, function(z) z %in% names(data) && anyNA(data[[z]][bad]), logical(1))]
    cl <- if(what == "q_r") data$cell[bad] else unique(c(data$from_cell[bad], data$cell[bad]))
    stop("Covariate(s) ", paste(v, collapse=", "), " in '", what, "' are NA for ", sum(bad), " row(s) ",
         "(cells ", paste(utils::head(sort(unique(cl)), 10), collapse=", "), if(length(unique(cl)) > 10) ", ...",
         "). Mask those cells with 'rast_mask' in make_walk_data() or fill the missing values.", call. = FALSE)
  }
  X
}

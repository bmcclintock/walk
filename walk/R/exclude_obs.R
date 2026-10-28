#' @title Exclude observations from walk data
#' @description Marks observations as excluded: their rows of `walk_data$L` are emptied, so
#' they carry no location information, but their times are kept. This is the same as the
#' `exclude` and `vmax` arguments of \code{\link{make_walk_data}} (which call this function), for
#' use on an existing `walk_ddl` object, e.g. after a fit has reported `fit$zero_obs`.
#' @param walk_data A `walk_ddl` object from \code{\link{make_walk_data}}.
#' @param obs Observations to exclude: rows of `walk_data$times` (= rows of `walk_data$L`, =
#' `walk_data$times$obs`), e.g. `fit$zero_obs`, or a logical vector with one value per row.
#' @param vmax Optional maximum plausible speed, in map units per `time_unit` of the data (e.g.
#' m/h), for a speed filter applied after `obs` are excluded (see 'Details').
#' @param reason Character label stored in `walk_data$times$excluded` for `obs`.
#' @details
#' **Excluded observations.** Because the movement model is a continuous-time Markov chain, an
#' excluded observation is exactly equivalent to deleting it and measuring the next time step
#' from the previous retained observation (\code{\link{fit_ctmc}} starts each segment at its first
#' retained observation, so this also holds for excluded first observations).
#' \code{\link{predict_ctmc}} still returns a location distribution at its time, from the
#' movement model only. `walk_data$times$excluded` gives the reason (`NA` for retained
#' observations):
#' * `"no_cells"`: the observation has zero likelihood in every cell of the model domain (e.g.
#'   all its error neighbourhood is masked by `rast_mask`). \code{\link{make_walk_data}} labels these
#'   automatically.
#' * `reason` (default `"user"`): given in `obs` / `make_walk_data(exclude = )`.
#' * `"speed"`: removed by the `vmax` filter.
#'
#' **Speed filter.** Within each segment, the distance between two observations is the gap
#' between the bounding boxes of their likelihood supports (0 if they overlap; so location error is
#' accounted for), less a tolerance of one cell. Iteratively, the fastest step between consecutive
#' retained observations is found; if it exceeds `vmax`, whichever of its two observations
#' leaves the slower speed between its neighbours when removed is excluded, until no step exceeds
#' `vmax`.
#' @return `walk_data` with emptied `L` rows and the character column `times$excluded`.
#' @export
exclude_obs <- function(walk_data, obs=NULL, vmax=NULL, reason="user"){
  if(!inherits(walk_data, "walk_ddl")) stop("'walk_data' must be a 'walk_ddl' object from make_walk_data().")
  if(!is.character(reason) || length(reason) != 1 || is.na(reason)) stop("'reason' must be a single character string.")
  N <- nrow(walk_data$L)
  if(is.null(walk_data$times$excluded)) walk_data$times$excluded <- NA_character_
  # rows without location information that are not labelled yet (e.g. edited by hand)
  walk_data <- exclude_mark(walk_data, which(Matrix::rowSums(walk_data$L) == 0), "no_cells")
  if(!is.null(obs)){
    if(is.logical(obs)){
      if(length(obs) != N || anyNA(obs)) stop("A logical 'obs' must have one non-missing value per row of 'walk_data$times'.")
      obs <- which(obs)
    }
    if(!is.numeric(obs) || any(is.na(obs) | obs < 1 | obs > N | obs != round(obs)))
      stop("'obs' must be row numbers of 'walk_data$times' (1..", N, ") or a logical vector.")
    walk_data <- exclude_mark(walk_data, as.integer(obs), reason)
  }
  if(!is.null(vmax)){
    if(!is.numeric(vmax) || length(vmax) != 1 || !is.finite(vmax) || vmax <= 0) stop("'vmax' must be a single positive number.")
    walk_data <- exclude_mark(walk_data, speed_filter(walk_data, vmax), "speed")
  }
  walk_data
}

# empty the L rows 'obs' and label them (existing labels are kept)
#' @noRd
exclude_mark <- function(walk_data, obs, reason){
  obs <- sort(unique(obs))
  if(!length(obs)) return(walk_data)
  walk_data$L <- exclude_rows(walk_data$L, obs)
  walk_data$times$excluded[obs] <- ifelse(is.na(walk_data$times$excluded[obs]), reason, walk_data$times$excluded[obs])
  walk_data
}

# empty rows of a sparse matrix (also drops stored zeros)
#' @noRd
exclude_rows <- function(L, rows){
  if(!length(rows)) return(L)
  if(!inherits(L, "dgCMatrix")) L <- as(L, "dgCMatrix")
  r <- drop_rows_csc(L@p, L@i, L@x, seq_len(nrow(L)) %in% rows)
  L@p <- r$p; L@i <- r$i; L@x <- r$x
  L
}

# rows to exclude with the maximum-speed filter (retained rows only), in the order removed
#' @noRd
speed_filter <- function(walk_data, vmax){
  a <- speed_filter_input(walk_data)
  L <- walk_data$L
  if(!inherits(L, "dgCMatrix")) L <- as(L, "dgCMatrix")
  bb <- row_support_bbox(L@p, L@i, L@x, nrow(L), as.numeric(walk_data$q_r$x), as.numeric(walk_data$q_r$y))
  speed_filter_cpp(bb, a$tt, as.integer(a$seg), a$keep, vmax, a$buffer)
}

# inputs shared by speed_filter() and the R reference implementation
#' @noRd
speed_filter_input <- function(walk_data){
  N <- nrow(walk_data$L)
  h <- walk_data$q_m$h
  seg <- walk_data$times$seg
  if(is.null(seg)) seg <- rep(1L, N)
  dt <- walk_data$times$dt
  dt[c(TRUE, seg[-1] != seg[-N])] <- 0
  list(buffer = if(length(h)) min(h) else 0, seg = seg, tt = cumsum(dt),     # differences within a segment are valid
       keep = is.na(walk_data$times$excluded))
}

# Plain-R version of speed_filter() (the algorithm documented in exclude_obs()); used by the
# tests to check the compiled version.
#' @noRd
speed_filter_r <- function(walk_data, vmax){
  a <- speed_filter_input(walk_data)
  L <- walk_data$L
  N <- nrow(L)
  T <- Matrix::summary(as(L, "TsparseMatrix"))
  T <- T[T$x > 0, ]
  xs <- walk_data$q_r$x[T$j]; ys <- walk_data$q_r$y[T$j]
  bb <- matrix(NA_real_, N, 4)
  bb[, 1] <- tapply_full(xs, T$i, min, N); bb[, 2] <- tapply_full(xs, T$i, max, N)
  bb[, 3] <- tapply_full(ys, T$i, min, N); bb[, 4] <- tapply_full(ys, T$i, max, N)
  gap <- function(a, b){
    dx <- pmax(0, bb[a,1] - bb[b,2], bb[b,1] - bb[a,2])
    dy <- pmax(0, bb[a,3] - bb[b,4], bb[b,3] - bb[a,4])
    sqrt(dx^2 + dy^2)
  }
  tt <- a$tt; seg <- a$seg; keep <- a$keep & !is.na(bb[, 1]); buffer <- a$buffer   # rows without cells never take part
  flagged <- integer(0)
  spd <- function(a, b){                       # minimum speed needed between observations a and b
    d <- pmax(0, gap(a, b) - buffer); t <- abs(tt[b] - tt[a])
    ifelse(d == 0, 0, ifelse(t > 0, d / t, Inf))
  }
  for(s in unique(seg)){
    repeat{
      idx <- which(seg == s & keep)
      n <- length(idx)
      if(n < 2) break
      sp <- spd(idx[-n], idx[-1])                 # speed of step j -> j+1
      j <- which.max(sp)
      if(!(sp[j] > vmax)) break
      # drop whichever of the two observations leaves the slower bypass
      by_a <- if(j > 1) spd(idx[j - 1], idx[j + 1]) else 0          # remove idx[j]
      by_b <- if(j + 1 < n) spd(idx[j], idx[j + 2]) else 0          # remove idx[j+1]
      drop <- if(by_b < by_a) idx[j + 1] else if(by_a < by_b) idx[j] else
        idx[j + if(j + 1 == n) 1 else 0]          # tie: the earlier row, except on the last step (the last row)
      keep[drop] <- FALSE
      flagged <- c(flagged, drop)
    }
  }
  flagged
}

#' @noRd
tapply_full <- function(x, i, f, N){
  out <- rep(NA_real_, N)
  if(length(x)){ v <- tapply(x, i, f); out[as.integer(names(v))] <- v }
  out
}

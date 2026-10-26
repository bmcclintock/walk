# Small simulated data sets for the joint-likelihood tests.
# Everything is simulated; no external data or langevinSSM needed.
# Kept small (~1,200 cells, < 40 fixes per track) so the whole suite runs in
# a few seconds on CRAN, and never uses more than 2 threads.

sim_covs <- function(res = 2000, nx = 60e3, ny = 40e3, seed = 1) {
  set.seed(seed)
  r <- terra::rast(xmin = 0, xmax = nx, ymin = 0, ymax = ny, resolution = res,
                   crs = "EPSG:3571")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  c1 <- terra::rast(r); terra::values(c1) <- sin(xy[, 1] / 9000) + cos(xy[, 2] / 7000)
  c2 <- terra::rast(r); terra::values(c2) <- as.vector(scale(xy[, 1] + 0.5 * xy[, 2]))
  out <- c(c1, c2); names(out) <- c("c1", "c2")
  out
}

sim_track <- function(id, n, seed, start, sigma = 1500, box = c(0, 60e3, 0, 40e3),
                      t0 = as.POSIXct("2016-09-16", tz = "UTC")) {
  set.seed(seed)
  dt <- c(0, stats::rexp(n - 1, 1 / 1.5))            # hours
  x <- start[1]; y <- start[2]; X <- Y <- numeric(n)
  for (i in seq_len(n)) {
    if (i > 1) {
      x <- x + stats::rnorm(1, 0, sigma * sqrt(dt[i]))
      y <- y + stats::rnorm(1, 0, sigma * sqrt(dt[i]))
    }
    x <- min(max(x, box[1] + 8e3), box[2] - 8e3)
    y <- min(max(y, box[3] + 8e3), box[4] - 8e3)
    X[i] <- x; Y[i] <- y
  }
  gps  <- stats::runif(n) < 0.3
  smaj <- ifelse(gps, NA, pmin(stats::rlnorm(n, log(2500), 0.5), 5000))
  smin <- ifelse(gps, NA, smaj * stats::runif(n, 0.1, 0.4))
  eor  <- ifelse(gps, NA, stats::runif(n, 0, pi))
  d <- data.frame(id = id, date = t0 + 3600 * cumsum(dt),
                  x = X + stats::rnorm(n, 0, 300), y = Y + stats::rnorm(n, 0, 300),
                  smaj = smaj, smin = smin, eor = eor,
                  x.err = ifelse(gps, 50, NA), y.err = ifelse(gps, 50, NA))
  class(d) <- c("dataLangevin", "data.frame")
  d
}

as_dl <- function(d) { class(d) <- c("dataLangevin", "data.frame"); d }

mp_test <- function() {
  ctmc_control(q_r = ctmc_model(~ 1, link = "log"),
               q_m = ctmc_model(~ d_c1 + d_c2, link = "log"),
               norm = FALSE, delta = "uniform")
}

add_d <- function(wd) {
  wd$q_m$d_c1 <- wd$q_m$c1 - wd$q_m$from_c1
  wd$q_m$d_c2 <- wd$q_m$c2 - wd$q_m$from_c2
  wd
}

dl_of <- function(wd, ...) fit_ctmc(add_d(wd), model_parameters = mp_test(), fit = FALSE, ...)$data_list

par_test <- function() c(log(1500^2 / (2 * 2000^2)), 0.3, -0.2)

sim_world <- function() {
  cv <- sim_covs()
  A <- sim_track("A", 30, 11, c(25e3, 20e3))
  B <- sim_track("B", 25, 12, c(35e3, 22e3))
  B$date <- B$date - 30 * 86400          # B starts earlier: stacked time is not monotone
  list(cv = cv, A = A, B = B)
}

# Minimal stand-in for ctmm telemetry objects (ctmm is not installed in this sandbox).
# ctmm/R/1.R:  new.telemetry <- methods::setClass("telemetry", contains="data.frame",
#                representation=methods::representation(info="list", UERE="UERE"))
# Identity is stored in the 'info' slot (ctmm/R/telemetry.R: info <- list(identity=id[i], ...)).
# We mirror that: an S4 class extending data.frame with an 'info' slot, and the columns
# proc_telem() reads (timestamp, longitude, latitude, COV.x.x, COV.x.y, COV.y.y).
if(!methods::isClass("telemetry")){
  methods::setClass("telemetry", contains = "data.frame",
                    representation = methods::representation(info = "list"))
}

# Convert a simulated dataLangevin track (x/y in the CRS of 'cv') into a mock telemetry object.
# Error covariance: GPS fixes keep x.err/y.err; Argos fixes use the ellipse covariance.
mock_telemetry <- function(d, cv, identity = as.character(d$id[1])){
  cov <- argos_diag_to_cov(d$smaj, d$smin, 180 * d$eor / pi)
  g <- !is.na(d$x.err)
  cov$cov.x.x[g] <- d$x.err[g]^2; cov$cov.y.y[g] <- d$y.err[g]^2; cov$cov.x.y[g] <- 0
  ll <- terra::crds(terra::project(terra::vect(cbind(d$x, d$y), crs = terra::crs(cv)), "epsg:4326"))
  df <- data.frame(timestamp = d$date, t = as.numeric(d$date),
                   longitude = ll[, 1], latitude = ll[, 2],
                   COV.x.x = cov$cov.x.x, COV.x.y = cov$cov.x.y, COV.y.y = cov$cov.y.y)
  methods::new("telemetry", df, info = list(identity = identity, timezone = "UTC"))
}

# time_unit = NULL (default): use the time unit stored on the data object.
# Simulated data only (helper-sim.R).

with_unit <- function(d, u) { attr(d, "time.unit") <- u; d }
num_dates <- function(d, per = 3600) { d$date <- as.numeric(d$date - d$date[1]) / per; d }   # numeric dates

test_that("proc_Langevin uses attr(data, 'time.unit') by default and time_unit overrides it", {
  w <- sim_world()
  ds <- c(0, diff(as.numeric(w$A$date)))                          # seconds
  # no attribute, POSIXct dates -> hours (formatData()'s default time.unit)
  p0 <- suppressWarnings(proc_Langevin(w$A, w$cv))
  expect_equal(attr(p0$times, "time_unit"), "hours"); expect_equal(p0$times$dt, ds / 3600)
  # attribute set by langevinSSM::formatData(time.unit = ...)
  for (u in c("secs", "mins", "hours", "days")) {
    p <- suppressWarnings(proc_Langevin(with_unit(w$A, u), w$cv))
    m <- c(secs = 1, mins = 60, hours = 3600, days = 86400)[[u]]
    expect_equal(attr(p$times, "time_unit"), u); expect_equal(p$times$dt, ds / m)
    expect_equal(p$L, p0$L)
  }
  # explicit time_unit overrides the attribute (and converts)
  p <- suppressWarnings(proc_Langevin(with_unit(w$A, "mins"), w$cv, time_unit = "days"))
  expect_equal(attr(p$times, "time_unit"), "days"); expect_equal(p$times$dt, ds / 86400)
  # aliases
  p <- suppressWarnings(proc_Langevin(with_unit(w$A, "minutes"), w$cv))
  expect_equal(attr(p$times, "time_unit"), "mins")
  expect_error(suppressWarnings(proc_Langevin(with_unit(w$A, "fortnights"), w$cv)), "time unit")
  # max_gap is in the resolved unit
  pg <- suppressWarnings(proc_Langevin(with_unit(w$A, "mins"), w$cv, max_gap = 120))
  pg2 <- suppressWarnings(proc_Langevin(w$A, w$cv, max_gap = 2))
  expect_identical(pg$times$seg, pg2$times$seg)
})

test_that("numeric dates (simLangevin) are used as they are", {
  w <- sim_world()
  A <- num_dates(w$A)                       # numeric, in hours, no attribute (like simLangevin)
  dh <- c(0, diff(A$date))
  p <- suppressWarnings(proc_Langevin(A, w$cv))
  expect_equal(attr(p$times, "time_unit"), "native"); expect_equal(p$times$dt, dh)
  # time_unit names the unit of the numeric dates, no conversion
  p <- suppressWarnings(proc_Langevin(A, w$cv, time_unit = "hours"))
  expect_equal(attr(p$times, "time_unit"), "hours"); expect_equal(p$times$dt, dh)
  # numeric dates with an attribute are in that unit and can be converted
  p <- suppressWarnings(proc_Langevin(with_unit(A, "hours"), w$cv))
  expect_equal(attr(p$times, "time_unit"), "hours"); expect_equal(p$times$dt, dh)
  p <- suppressWarnings(proc_Langevin(with_unit(A, "hours"), w$cv, time_unit = "mins"))
  expect_equal(attr(p$times, "time_unit"), "mins"); expect_equal(p$times$dt, dh * 60)
  # same likelihood as the POSIXct version in hours
  wdA <- make_walk_data(suppressWarnings(proc_Langevin(A, w$cv)), w$cv)
  wdP <- make_walk_data(suppressWarnings(proc_Langevin(w$A, w$cv)), w$cv)
  par <- par_test()
  expect_equal(as.numeric(ctmc_n2ll_joint(par, dl_of(wdA))), as.numeric(ctmc_n2ll_joint(par, dl_of(wdP))),
               tolerance = 1e-12)
  # several animals with numeric dates: dt within animal, no negative steps
  B <- num_dates(w$B)
  pAB <- suppressWarnings(proc_Langevin(as_dl(rbind(A, B)), w$cv))
  expect_true(all(pAB$times$dt >= 0)); expect_equal(unique(pAB$times$seg), 1:2)
  expect_equal(pAB$times$dt, c(dh, 0, diff(B$date)))
  # mixed inputs are rejected
  expect_error(suppressWarnings(proc_Langevin(list(A, w$B), w$cv)), "different classes")
  expect_error(suppressWarnings(proc_Langevin(list(with_unit(w$A, "mins"), with_unit(w$B, "hours")), w$cv)),
               "different time units")
})

test_that("proc_telem: time_unit defaults to hours, carries ts_scale", {
  w <- sim_world()
  tA <- mock_telemetry(w$A, w$cv)
  ds <- c(0, diff(as.numeric(w$A$date)))
  p <- suppressWarnings(proc_telem(tA, w$cv))
  expect_equal(attr(p$times, "time_unit"), "hours"); expect_equal(p$times$dt, ds / 3600)
  expect_equal(attr(p$times, "ts_scale"), 1 / 3600)
  expect_identical(p$times, suppressWarnings(proc_Langevin(w$A, w$cv))$times)
  for (u in c("secs", "mins", "days")) {
    m <- c(secs = 1, mins = 60, days = 86400)[[u]]
    p <- suppressWarnings(proc_telem(tA, w$cv, time_unit = u))
    expect_equal(attr(p$times, "time_unit"), u); expect_equal(p$times$dt, ds / m)
    expect_equal(attr(p$times, "ts_scale"), 1 / m)
  }
  expect_error(suppressWarnings(proc_telem(tA, w$cv, time_unit = NULL)), "time_unit")
  # max_gap is in time_unit
  expect_identical(suppressWarnings(proc_telem(tA, w$cv, time_unit = "secs", max_gap = 7200))$times$seg,
                   suppressWarnings(proc_telem(tA, w$cv, max_gap = 2))$times$seg)
})

test_that("bind_proc_data and predict_ctmc carry the time scale", {
  w <- sim_world()
  pA <- suppressWarnings(proc_Langevin(with_unit(w$A, "mins"), w$cv))
  pB <- suppressWarnings(proc_telem(mock_telemetry(w$B, w$cv), w$cv, time_unit = "mins"))
  b <- bind_proc_data(pA, pB)
  expect_equal(attr(b$times, "time_unit"), "mins"); expect_equal(attr(b$times, "ts_scale"), 1 / 60)
  # POSIXct and numeric time stamps in the same nominal unit cannot be mixed
  pN <- suppressWarnings(proc_Langevin(num_dates(w$B, 60), w$cv, time_unit = "mins"))
  expect_error(bind_proc_data(pA, pN), "different scales")
  # predict_ctmc: dt at the observation times equals the proc dt (POSIXct and numeric dates)
  for (A in list(with_unit(w$A, "mins"), num_dates(w$A))) {
    wd <- add_d(make_walk_data(suppressWarnings(proc_Langevin(A, w$cv)), w$cv))
    f <- fit_ctmc(wd, model_parameters = mp_test(), fit = FALSE, start = list(beta_q_r = par_test()[1], beta_q_m = par_test()[-1]))
    expect_equal(f$time_unit, attr(wd$times, "time_unit"))
    pr <- predict_ctmc(f, wd)
    o <- pr$times[pr$times$type == "o", ]
    expect_equal(o$dt, wd$times$dt, tolerance = 1e-12)
    expect_equal(attr(pr$times, "time_unit"), attr(wd$times, "time_unit"))
  }
})

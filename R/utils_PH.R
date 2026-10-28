# Scenario 1 (proportional hazards): data generation, amputation, and Cox estimation.
# h(t | x) = lambda * exp(0.8 x1 + 0.5 x2 + 0.2 x3); x2 and x3 are made incomplete.

library(survival)
library(mice)

#' True log hazard ratios of x1, x2, x3 in Scenario 1
TRUE_LOGHR_PH <- c(x1 = 0.8, x2 = 0.5, x3 = 0.2)

#' Generate one complete dataset of Scenario 1
#'
#' Event times are drawn piecewise on a fine time grid (step 0.001) with a constant
#' hazard ratio; censoring is exponential, with administrative censoring at maxt.
#'
#' @param n sample size
#' @param maxt administrative censoring time (years)
#' @param hr true log hazard ratios of x1, x2, x3
#' @return data.frame with time, status, x1, x2, x3 and the Nelson-Aalen cumulative
#'   hazard (cumhaz, used as an auxiliary variable by the imputation models)
genPH <- function(n = 1000, maxt = 5, hr = TRUE_LOGHR_PH) {
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5)
  lambda_c <- 0.08   # censoring hazard, about 30% censored
  lambda <- 0.2      # baseline hazard
  delta <- 0.001
  t_grid <- seq(0, maxt - delta, delta)
  t_event_pw <- matrix(nrow = n, ncol = length(t_grid))
  for (i in seq_along(t_grid)) {
    u <- runif(n)
    t_a <- t_grid[i] - log(u) / (lambda * exp(hr[1] * x1 + hr[2] * x2 + hr[3] * x3))
    t_event_pw[, i] <- ifelse(t_a >= t_grid[i] & t_a < t_grid[i] + delta, t_a, NA)
  }
  t_event <- suppressWarnings(apply(t_event_pw, 1, min, na.rm = TRUE))
  t_event[is.infinite(t_event)] <- 10
  t_cens <- -log(runif(n)) / lambda_c
  status <- as.integer(t_event < t_cens & pmin(t_event, t_cens) <= maxt)
  time <- pmin(t_event, t_cens, maxt)
  dat <- data.frame(time, status, x1, x2, x3)
  dat$cumhaz <- mice::nelsonaalen(dat, time, status)
  dat
}

#' Introduce missing values in x2 and x3 with mice::ampute()
#'
#' @param dat complete data from genPH()
#' @param mech "MCAR", "MAR" or "MNAR"
#' @param pmiss proportion of incomplete rows
#' @return amputed data.frame, with x3 as a factor
generate_NA_PH <- function(dat, mech, pmiss = 0.3) {
  stopifnot(mech %in% c("MCAR", "MAR", "MNAR"))
  patterns <- as.integer(!names(dat) %in% c("x2", "x3"))
  amp <- suppressWarnings(mice::ampute(dat, mech = mech, prop = pmiss,
                                       patterns = patterns)$amp)
  amp$x3 <- factor(amp$x3)
  amp
}

#' Performance of one Cox fit: estimate, SE, 95% CI, bias, relative bias, squared
#' error and coverage of the true log hazard ratio, for each covariate
#'
#' @param est,se estimates and standard errors of x1, x2, x3
#' @param crit critical value(s) of the 95% interval (normal, or t for Rubin's rules)
#' @param truelogHR true log hazard ratios
perf_PH <- function(est, se, crit = qnorm(0.975), truelogHR = TRUE_LOGHR_PH) {
  lo95 <- est - crit * se
  hi95 <- est + crit * se
  data.frame(covariates = names(truelogHR), est = est, se_est = se, lo95 = lo95,
             hi95 = hi95, bias = est - truelogHR,
             bias_rel = (est - truelogHR) / truelogHR * 100,
             MSE_estimate = (est - truelogHR)^2,
             cover = truelogHR >= lo95 & truelogHR <= hi95, row.names = NULL)
}

#' Cox model on one complete or singly imputed dataset
#'
#' @param data data.frame with time, status, x1, x2, x3
evaluate_coxest_PH <- function(data) {
  fit <- coxph(Surv(time, status) ~ x1 + x2 + x3, data = data)
  perf_PH(unname(coef(fit)), unname(sqrt(diag(vcov(fit)))))
}

#' Cox model on each multiply imputed dataset, pooled with Rubin's rules
#' (t reference distribution with the Barnard-Rubin degrees of freedom)
#'
#' @param mids mids object returned by mice
evaluate_coxest_PH_multiple <- function(mids) {
  fit <- with(mids, coxph(Surv(time, status) ~ x1 + x2 + x3))
  res <- summary(mice::pool(fit))
  perf_PH(res$estimate, res$std.error, crit = qt(0.975, res$df))
}

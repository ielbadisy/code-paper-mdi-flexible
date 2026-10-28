# Scenarios 2-4 (time-varying effect of x3): data generation, amputation, estimation.
# h(t | x) = lambda * exp(0.8 x1 + 0.5 x2 + beta(t) x3), with beta(t) given by the shape.

library(survival)
library(rstpm2)
library(mice)

#' True time-varying log hazard ratio beta(t) of x3 for each shape (Table 1)
TVE_BETA <- list(
  bathub = function(t) (0.01 + 1.01 * exp(-t / 0.3) + 0.06 * (t^1.1)),  # Scenario 2
  decreasing = function(t) (1.2 / exp(t) - 0.02 * t^0.7),               # Scenario 3
  increasing = function(t) (0.01 + 0.65 * t^0.4)                        # Scenario 4
)

#' Constant added to the linear predictor (not multiplied by x3), by shape
TVE_OFFSET <- c(bathub = 0, decreasing = 0.01, increasing = 0)

#' Degrees of freedom of the time spline of x3 in the estimation model, by shape
TVE_DF <- c(bathub = 4, decreasing = 2, increasing = 2)

#' True hazard ratio of x3 at time t
#'
#' @param t time (years)
#' @param shape "bathub", "decreasing" or "increasing"
true_hr_TVE <- function(t, shape) exp(TVE_BETA[[shape]](t))

#' Generate one complete dataset of Scenarios 2-4
#'
#' Event times are drawn piecewise on a time grid (step 0.01), with the hazard
#' ratio of x3 evaluated at the start of each interval; censoring is exponential,
#' with administrative censoring at maxt.
#'
#' @param n sample size
#' @param maxt administrative censoring time (years)
#' @param shape "bathub", "decreasing" or "increasing"
#' @return data.frame with time, status, x1, x2, x3 and the Nelson-Aalen cumulative
#'   hazard (cumhaz, auxiliary variable for imputation)
genTVE <- function(n = 1000, maxt = 5, shape) {
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5)
  lambda_c <- 0.15   # censoring hazard, about 30% censored
  lambda <- 0.2      # baseline hazard
  beta <- TVE_BETA[[shape]]
  lp <- 0.8 * x1 + 0.5 * x2
  if (TVE_OFFSET[[shape]] != 0) lp <- lp + TVE_OFFSET[[shape]]
  delta <- 0.01
  t_grid <- seq(0, maxt - delta, delta)
  t_event_pw <- matrix(nrow = n, ncol = length(t_grid))
  for (i in seq_along(t_grid)) {
    u <- runif(n)
    t_a <- t_grid[i] - log(u) / (lambda * exp(lp + beta(t_grid[i]) * x3))
    t_event_pw[, i] <- ifelse(t_a >= t_grid[i] & t_a < t_grid[i] + delta, t_a, NA)
  }
  t_event <- suppressWarnings(apply(t_event_pw, 1, min, na.rm = TRUE))
  t_event[is.infinite(t_event)] <- 30
  t_cens <- -log(runif(n)) / lambda_c
  status <- as.integer(t_event < t_cens & pmin(t_event, t_cens) <= maxt)
  time <- pmin(t_event, t_cens, maxt)
  dat <- data.frame(time, status, x1, x2, x3)
  cumhaz <- mice::nelsonaalen(dat, time, status)
  dat$cumhaz <- ifelse(is.na(cumhaz), 0, cumhaz)
  dat
}

#' Introduce missing values in x3 with mice::ampute()
#'
#' @param dat complete data from genTVE()
#' @param mech "MCAR", "MAR" or "MNAR"
#' @param pmiss proportion of incomplete rows
generate_NA_TVE <- function(dat, mech, pmiss = 0.3) {
  stopifnot(mech %in% c("MCAR", "MAR", "MNAR"))
  patterns <- as.integer(names(dat) != "x3")
  suppressWarnings(mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp)
}

#' Fit the flexible parametric model and return the hazard ratio curve of x3
#'
#' The curve is predicted on the default rstpm2 time grid (300 points over the
#' observed follow-up), with 95% confidence limits.
#'
#' @param dat complete or imputed data
#' @param df degrees of freedom of the time spline of x3
#' @return matrix with columns time, hr, lower, upper
fit_TVE <- function(dat, df) {
  fit <- stpm2(Surv(time, status) ~ x1 + x2 + x3, data = dat, tvc = list(x3 = df))
  p <- predict(fit, type = "hr", newdata = data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",
               grid = TRUE, full = TRUE, se.fit = TRUE)
  cbind(time = p$time, hr = p$Estimate, lower = p$lower, upper = p$upper)
}

#' Pool the hazard ratio curves of m imputed datasets with Rubin's rules
#'
#' Pointwise on the hazard ratio scale: the pooled curve is the mean of the m
#' curves; the within-imputation variance is derived from the width of each 95%
#' interval; the 95% interval uses the total variance with a normal reference.
#'
#' @param fits list of matrices returned by fit_TVE()
#' @return matrix with columns time, hr, lower, upper
pool_TVE <- function(fits) {
  a <- simplify2array(fits)
  hr <- a[, "hr", ]
  within <- rowMeans(((a[, "upper", ] - a[, "lower", ]) / (2 * 1.96))^2)
  total <- within + (1 + 1 / length(fits)) * apply(hr, 1, var)
  q <- rowMeans(hr)
  cbind(time = rowMeans(a[, "time", ]), hr = q,
        lower = q - 1.96 * sqrt(total), upper = q + 1.96 * sqrt(total))
}

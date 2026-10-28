#====================[ SCENARIO 5: NON-LINEAR EFFECT OF x2 ]====================
# Data generation, amputation, estimation, and Rubin pooling for Scenario 5.
# The effect of x2 is non-linear and constant over time:
#   h(t | x) = lambda * exp(0.8 x1 + g(x2) + 0.3 x3),  g(x2) = 0.6 log(1 + (x2 - 1.5)^2)
# The estimand is the log hazard ratio g(x) - g(1.5) on a grid of x2 values.

# one thread per worker: parallelism comes from replicates, not from the imputers
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(ranger.num.threads = 1)

suppressPackageStartupMessages({
  library(survival); library(rstpm2); library(mice); library(splines)
})

LAMBDA   <- 0.2    # baseline hazard
# censoring hazard: about 30% censored, including administrative censoring at MAXT
LAMBDA_C <- 0.058
MAXT     <- 5      # administrative censoring (years)

g_x2  <- function(x) 0.6 * log(1 + (x - 1.5)^2)   # true log hazard ratio, minimum at x2 = 1.5
X_REF <- 1.5                                       # reference value of x2

#---------- spline basis with knots fixed from the population distribution of x2 ----------
NS_DF <- 7
local({
  set.seed(1)
  z  <- rnorm(1e6)
  xp <- z^2 + z + runif(1e6)
  NS_KNOTS <<- quantile(xp, seq_len(NS_DF - 1) / NS_DF)
  NS_BOUND <<- quantile(xp, c(0.025, 0.975))
  X_GRID   <<- seq(quantile(xp, 0.05), quantile(xp, 0.95), length.out = 25)
})
BN  <- paste0("b", seq_len(NS_DF))
nsb <- function(x) {
  B <- ns(x, knots = NS_KNOTS, Boundary.knots = NS_BOUND)
  colnames(B) <- BN
  B
}

#---------- generate complete data ----------
# Nelson-Aalen cumulative hazard at each observed time. Same estimator as
# mice::nelsonaalen(), which can return NA through a floating-point time mismatch.
nelson_aalen <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  stats::stepfun(fit$time, c(0, cumsum(fit$n.event / fit$n.risk)))(time)
}

genNLE <- function(n = 1000) {
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5)
  t_event <- rexp(n, LAMBDA * exp(0.8 * x1 + g_x2(x2) + 0.3 * x3))
  t_cens  <- rexp(n, LAMBDA_C)
  time    <- pmin(t_event, t_cens, MAXT)
  status  <- as.integer(t_event <= pmin(t_cens, MAXT))
  dat <- data.frame(time, status, x1, x2, x3)
  dat$cumhaz <- nelson_aalen(time, status)   # auxiliary variable for imputation
  dat
}

#---------- introduce NA in x2 only ----------
generate_NA_NLE <- function(dat, mech, pmiss = 0.3) {
  patterns <- as.integer(names(dat) != "x2")
  suppressWarnings(mice::ampute(dat, prop = pmiss, mech = mech, patterns = patterns)$amp)
}

#---------- fit the flexible model and return log-HR and SE on the x2 grid ----------
fit_NLE <- function(d) {
  d   <- cbind(d, nsb(d$x2))
  fit <- stpm2(reformulate(c("x1", BN, "x3"), "Surv(time, status)"), data = d)
  b   <- coef(fit)[BN]
  V   <- vcov(fit)[BN, BN]
  C   <- sweep(nsb(X_GRID), 2, nsb(X_REF))               # contrast with x2 = 1.5
  data.frame(x2 = X_GRID, est = drop(C %*% b), se = sqrt(rowSums((C %*% V) * C)))
}

#---------- confidence intervals ----------
ci_single <- function(f) {
  data.frame(x2 = f$x2, est = f$est,
             lower = f$est - qnorm(0.975) * f$se, upper = f$est + qnorm(0.975) * f$se)
}

ci_rubin <- function(fits) {           # Rubin's rules on the log-HR scale, per grid point
  E <- sapply(fits, `[[`, "est"); S <- sapply(fits, `[[`, "se"); m <- ncol(E)
  Q <- rowMeans(E); W <- rowMeans(S^2); B <- apply(E, 1, var)
  Tv <- W + (1 + 1 / m) * B
  nu <- (m - 1) * (1 + W / ((1 + 1 / m) * B))^2
  q  <- qt(0.975, nu)
  data.frame(x2 = fits[[1]]$x2, est = Q, lower = Q - q * sqrt(Tv), upper = Q + q * sqrt(Tv))
}

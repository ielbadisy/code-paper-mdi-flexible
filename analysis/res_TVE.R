# Scenarios 2-4: Figures 1-3 of the paper and Figures B1-B24 of Appendix B.
# Input: results/sim/res_TVE_<shape>.rds from sim/sim_TVE.R.
# Output: results/figures/Figure{1,2,3}_*.png, results/figures/appendix_B/*.png,
#         results/tables/TVE_metrics.csv (coverage and MSE at each evaluation time).

source("R/labels.R")
source("R/utils_TVE.R")
library(data.table)
library(ggplot2)

EVAL_TIMES <- seq(0.5, 5, by = 0.5)   # evaluation times of coverage and MSE
CURVE_TIMES <- seq(0.05, 5, by = 0.05) # common grid used to average the curves
SHAPES <- c("bathub", "decreasing", "increasing")
MECHS <- c("MCAR", "MAR", "MNAR")

#' Interpolate each replicate's curve at given times
#'
#' Every fitted curve is predicted on its own time grid (the observed follow-up of the
#' replicate), so curves are linearly interpolated at common times before averaging.
#' @param d results of one shape (data.table)
#' @param times times at which to interpolate
interp_curves <- function(d, times) {
  d[, {
    f <- function(y) approx(time, y, xout = times, rule = 2)$y
    list(time = times, hr = f(hr), lower = f(lower), upper = f(upper))
  }, by = .(mech, method, rep)]
}

#' Coverage(t) and MSE(t) of the hazard ratio of x3 against the true hazard ratio
tve_metrics <- function(pts, shape) {
  pts[, true_hr := true_hr_TVE(time, shape)]
  pts[, .(coverage = mean(lower <= true_hr & true_hr <= upper),
          mse = mean((hr - true_hr)^2), n_rep = .N),
      by = .(mech, method, time)]
}

plot_hr <- function(d, avg, shape, mech_) {
  ggplot(d[mech == mech_], aes(time, hr)) +
    geom_line(aes(group = rep), colour = "grey", linewidth = 0.2) +
    geom_function(fun = function(t) true_hr_TVE(t, shape), colour = "black") +
    geom_line(data = avg[mech == mech_], colour = "red", linetype = "11", linewidth = 0.8) +
    facet_wrap(vars(label)) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(x = "Time (years)", y = expression(HR[X[3]])) +
    theme_bw()
}

plot_coverage <- function(met, mech_) {
  ggplot(met[mech == mech_], aes(time, coverage)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_point(shape = 17, size = 2, colour = "red") +
    facet_wrap(vars(label), ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:5) +
    labs(x = "Time (years)", y = expression("Coverage rate of HR"[X[3]])) +
    theme_bw()
}

plot_mse <- function(pts, met, mech_) {
  ggplot(pts[mech == mech_], aes(time, (hr - true_hr)^2)) +
    geom_boxplot(aes(group = time), fill = "lightblue", outlier.colour = "grey40",
                 outlier.size = 0.1) +
    geom_point(data = met[mech == mech_], aes(y = mse), shape = 17, size = 1.5,
               colour = "red") +
    facet_wrap(vars(label), ncol = 4) +
    coord_cartesian(ylim = c(0, 3)) +
    scale_x_continuous(breaks = 1:5) +
    labs(x = "Time (years)", y = expression(paste("Squared error of ", HR[X[3]]))) +
    theme_bw()
}

dir.create(file.path("results", "figures", "appendix_B"), recursive = TRUE,
           showWarnings = FALSE)
all_metrics <- list()
for (shape in SHAPES) {
  d <- as.data.table(readRDS(file.path("results", "sim", sprintf("res_TVE_%s.rds", shape))))
  d <- d[method != "complete_data"]
  d[, label := method_factor(method)]
  avg <- interp_curves(d, CURVE_TIMES)[, .(hr = mean(hr)), by = .(mech, method, time)]
  avg[, label := method_factor(method)]
  pts <- interp_curves(d, EVAL_TIMES)
  met <- tve_metrics(pts, shape)
  pts[, label := method_factor(method)]
  met[, label := method_factor(method)]
  all_metrics[[shape]] <- cbind(shape = shape, met[, !"label"])
  for (m in MECHS) {
    out <- c(sprintf("TVE_%s_%s.png", m, shape), sprintf("coverage_TVE_%s_%s.png", shape, m),
             sprintf("MSE_TVE_%s_%s_boxplot.png", shape, m))
    plots <- list(plot_hr(d, avg, shape, m), plot_coverage(met, m), plot_mse(pts, met, m))
    for (k in 1:3) save_figure(plots[[k]], file.path("appendix_B", out[k]))
    if (shape == "increasing" && m == "MNAR") {
      save_figure(plots[[1]], "Figure1_TVE_increasing_MNAR.png")
      save_figure(plots[[2]], "Figure2_coverage_TVE_increasing_MNAR.png")
      save_figure(plots[[3]], "Figure3_MSE_TVE_increasing_MNAR.png")
    }
  }
  rm(d, pts); gc()
}
fwrite(rbindlist(all_metrics), results_path("tables", "TVE_metrics.csv"))

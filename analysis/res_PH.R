# Scenario 1 (proportional hazards): performance tables and Appendix A figures.
# Reads results/sim/res_PH.rds (from sim/sim_PH.R) and writes
#   results/tables/PH_performance.csv  (Tables 2-3: MNAR rows; Appendix A tables: all rows)
#   results/tables/table_2.tex, table_3.tex
#   results/figures/appendix_A/PH_bias.png, PH_Emp.png, PH_mse.png
# Run from the repository root: Rscript analysis/res_PH.R [path to res_PH.rds]
suppressPackageStartupMessages({ library(ggplot2) })
source("R/labels.R")

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else results_path("sim", "res_PH.rds")
res <- if (grepl("\\.csv$", path)) read.csv(path) else readRDS(path)

# x1 is a fully observed covariate; the tables and figures report x2 and x3
COVARIATES <- c(x2 = "X2 continuous", x3 = "X3 binary", x31 = "X3 binary")
res <- res[res$covariates %in% names(COVARIATES), ]
res$covariates <- factor(COVARIATES[res$covariates], levels = unique(COVARIATES))
res$method <- method_factor(res$method)
res$mech <- factor(res$mech, levels = c("MCAR", "MAR", "MNAR"))

#---------- performance averaged over replicates ----------
perf <- aggregate(cbind(bias, bias_rel, MSE_estimate, cover, se_est) ~
                    covariates + mech + method, data = res, FUN = mean)
names(perf) <- c("covariate", "mechanism", "method", "bias", "relative_bias_pct", "MSE",
                 "coverage", "mean_SE")
perf$replicates <- aggregate(bias ~ covariates + mech + method, data = res, FUN = length)$bias
perf <- perf[order(perf$covariate, perf$mechanism, perf$method), ]
write.csv(perf, results_path("tables", "PH_performance.csv"), row.names = FALSE)

#---------- Tables 2 and 3 (MNAR) ----------
f3 <- function(x) formatC(x, format = "f", digits = 3)
write_table <- function(cov, file) {
  d <- perf[perf$covariate == cov & perf$mechanism == "MNAR", ]
  lab <- sub("CompleteCases", "Complete cases", as.character(d$method))
  writeLines(c("Method & Bias & Relative bias (\\%) & MSE & Empirical coverage \\\\", "\\midrule",
               sprintf("%s & %s & %s & %s & %s \\\\", lab, f3(d$bias),
                       f3(d$relative_bias_pct), f3(d$MSE), f3(d$coverage))),
             results_path("tables", file))
}
write_table("X2 continuous", "table_2.tex")
write_table("X3 binary", "table_3.tex")

#---------- Appendix A figures: replicate-level distributions ----------
dir.create(file.path("results", "figures", "appendix_A"), recursive = TRUE, showWarnings = FALSE)
box <- function(var, xlim = NULL, vline = TRUE) {
  p <- ggplot(res, aes(x = .data[[var]], group = method)) +
    geom_boxplot(aes(fill = method), outlier.shape = NA) +
    guides(fill = guide_legend(reverse = TRUE)) +
    facet_grid(mech ~ covariates) +
    theme_bw() +
    theme(axis.text.y = element_blank())
  if (vline) p <- p + geom_vline(xintercept = 0, linetype = 20, colour = "red")
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim, expand = FALSE)
  p
}
save_figure(box("bias"), file.path("appendix_A", "PH_bias.png"))
save_figure(box("se_est", xlim = c(0, 0.15), vline = FALSE), file.path("appendix_A", "PH_Emp.png"))
save_figure(box("MSE_estimate", xlim = c(0, 0.3)), file.path("appendix_A", "PH_mse.png"))

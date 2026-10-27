#====================[ SCENARIO 5: RESULTS ]====================
# Reads res_NLE.rds (from sim_NLE.R) and produces the effect, coverage, and MSE
# figures for each missingness mechanism, plus a performance table with Monte
# Carlo standard errors. Estimand: log HR g(x2) - g(1.5) on the x2 grid.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

res <- as.data.table(readRDS("res_NLE.rds"))
dir.create("figures", showWarnings = FALSE)

LAB <- c(complete_data = "FullData", complete = "CompleteCases", naive = "Median/mode",
         hotdeck = "Hotdeck", mice = "MICE", cart = "CART", famd = "FAMD", knn = "KNN",
         misscforest = "missCforest", missforest = "missForest", missranger = "missRanger",
         micecart = "miceCART", micerf = "miceRF")

#---------- failures (fits that returned NA) ----------
fail <- res[, .(failed_reps = uniqueN(ID[is.na(est)])), by = .(mech, method)]
res  <- res[!is.na(est)]
res[, method := factor(LAB[method], levels = LAB)]
res[, `:=`(cover = truth >= lower & truth <= upper, sqerr = (est - truth)^2)]

#---------- performance per method and mechanism, averaged over the x2 grid ----------
# MCSE of coverage: sqrt(p(1-p)/R); of bias: SD(est)/sqrt(R); of MSE: SD(sqerr)/sqrt(R).
pt <- res[, .(R = .N, bias = mean(est - truth), empse = sd(est), mse = mean(sqerr),
              cover = mean(cover),
              mcse_bias = sd(est) / sqrt(.N), mcse_mse = sd(sqerr) / sqrt(.N),
              mcse_cover = sqrt(mean(cover) * (1 - mean(cover)) / .N)),
          by = .(mech, method, x2)]
perf <- pt[, .(R = min(R), mean_abs_bias = mean(abs(bias)), max_abs_bias = max(abs(bias)),
               mean_mse = mean(mse), mean_cover = mean(cover), min_cover = min(cover),
               max_mcse_bias = max(mcse_bias), max_mcse_mse = max(mcse_mse),
               max_mcse_cover = max(mcse_cover)), by = .(mech, method)]
perf <- merge(perf, fail[, .(mech, method = factor(LAB[method], levels = LAB), failed_reps)],
              by = c("mech", "method"))
setorder(perf, mech, method)
fwrite(perf, "performance_NLE.csv")
fwrite(pt, "performance_NLE_by_x2.csv")
pm <- res[, .(prop_missing = mean(prop_missing)), by = mech]
print(pm)
print(perf[, .(mech, method, R, failed_reps, bias = round(mean_abs_bias, 3),
               mse = round(mean_mse, 4), cover = round(mean_cover, 3),
               min_cover = round(min_cover, 3), mcse_cover = round(max_mcse_cover, 3))])

#---------- figures (12 panels: the full-data reference is in the table only) ----------
th <- theme_bw(base_size = 14)
for (mk in c("MNAR", "MAR", "MCAR")) {
  d <- res[mech == mk & method != "FullData"]
  d[, method := droplevels(method)]
  m <- d[, .(est = mean(est), truth = truth[1]), by = .(method, x2)]
  p_eff <- ggplot(d, aes(x2, exp(est))) +
    geom_line(aes(group = ID), colour = "grey", linewidth = 0.3) +
    geom_line(data = m, aes(y = exp(truth)), colour = "black", linewidth = 0.7) +
    geom_line(data = m, colour = "red", linetype = "dotted", linewidth = 1.2) +
    facet_wrap(~method, ncol = 4) +
    labs(x = expression(X[2]), y = expression(HR[X[2]] ~ "(reference" ~ X[2] == 1.5 * ")")) + th
  p_cov <- ggplot(pt[mech == mk & method != "FullData"], aes(x2, cover)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_point(shape = 17, colour = "red", size = 2.5) +
    facet_wrap(~droplevels(method), ncol = 4) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = expression(X[2]), y = expression("Coverage rate of" ~ log ~ HR[X[2]])) + th
  # per-replicate MSE within each rounded x2 bin
  b <- d[, .(mse = mean(sqerr)), by = .(method, ID, bin = factor(round(x2)))]
  p_mse <- ggplot(b, aes(bin, mse)) +
    geom_boxplot(fill = "lightblue", outlier.colour = "red", outlier.size = 0.5) +
    facet_wrap(~method, ncol = 4) +
    labs(x = expression(X[2] ~ "(rounded)"), y = expression("MSE of" ~ log ~ HR[X[2]])) + th
  for (p in list(list("effect", p_eff), list("coverage", p_cov), list("mse", p_mse)))
    ggsave(sprintf("figures/NLE_%s_%s.png", p[[1]], mk), p[[2]], width = 12, height = 9, dpi = 200)
}

#---------- LaTeX table for Appendix C ----------
f3 <- function(x) formatC(x, format = "f", digits = 3)
tb <- perf[order(match(mech, c("MNAR", "MAR", "MCAR")), method)]
rows <- tb[, sprintf("%s & %s & %s & %s & %s & %s & %s & %d \\\\",
                     mech, sub("CompleteCases", "Complete cases", sub("FullData", "Full data", method)), f3(mean_abs_bias), f3(mean_mse),
                     f3(mean_cover), f3(min_cover), f3(max_mcse_cover), failed_reps)]
mid <- which(diff(match(tb$mech, c("MNAR", "MAR", "MCAR"))) != 0)
rows[mid] <- paste(rows[mid], "\\midrule")
writeLines(c(
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Performance in Scenario 5 averaged over the 25 grid values of $x_2$ (log hazard ratio scale, 500 replicates per mechanism). Mean $|$bias$|$: mean over the grid of the absolute bias. Min coverage: lowest coverage over the grid. MCSE: largest Monte Carlo standard error of coverage over the grid. Failed: replicates in which the method returned no estimate.}\\label{tab:c1}",
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "Mechanism & Method & Mean $|$bias$|$ & Mean MSE & Mean coverage & Min coverage & MCSE & Failed \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  "table_C1.tex")

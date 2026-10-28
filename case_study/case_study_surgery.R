# CRC application, Model 2: time-varying effect (TVE) of surgery (Figure 8).
# Run from the repository root: Rscript case_study/case_study_surgery.R
suppressMessages({
  library(survival); library(rstpm2); library(imputer); library(mice); library(miceadds)
  library(ggplot2)
})
source("R/labels.R")
source("R/utils_case_study.R")

CRC <- read_CRC()
DF_SURGERY <- 3  # df of the time-varying spline of surgery (see Appendix D)

#' Fit the TVE model and return HR(surgery) over follow-up time
evaluate_surgery <- function(dat) {
  fit <- stpm2(Surv(time, event) ~ Age + Surgery + MSI.MSS.status + Sex + Stage,
               data = dat, tvc = list(Surgery = DF_SURGERY))
  pr <- predict(fit, type = "hr", var = "Surgery",
                newdata = data.frame(Age = 60, Surgery = 0, MSI.MSS.status = 0,
                                     Sex = 0, Stage = 1),
                grid = TRUE, seqLength = 300, full = TRUE)
  data.frame(time = pr$time, hr = pr$Estimate)
}

set.seed(1234)
ampdat <- ampute_CRC(CRC, "Surgery")   # one amputed dataset, used by every method
reference <- evaluate_surgery(CRC)
res <- do.call(rbind, lapply(CASE_METHODS, impute_and_evaluate,
                             ampdat = ampdat, evaluate = evaluate_surgery))

# Mean absolute deviation from the reference curve on the HR scale (the HR is close
# to 0 early in follow-up, where the log scale would inflate small differences).
# Curves share the time grid because follow-up times are never amputed.
res$ref <- approx(reference$time, reference$hr, res$time, rule = 2)$y
mad <- aggregate(cbind(MAD = abs(hr - ref)) ~ method, data = res, FUN = mean)
mad$method <- method_factor(mad$method)
mad <- mad[order(mad$MAD), ]
print(mad, digits = 2, row.names = FALSE)

write.csv(rbind(data.frame(reference, method = "complete_data"), res[, names(res) != "ref"]),
          results_path("tables", "case_study_surgery_curves.csv"), row.names = FALSE)
write.csv(mad, results_path("tables", "case_study_surgery_MAD.csv"), row.names = FALSE)

res$method <- method_factor(res$method)
ref_panels <- merge(reference, data.frame(method = levels(res$method)))
ref_panels$method <- factor(ref_panels$method, levels = levels(res$method))
p <- ggplot() +
  geom_line(data = ref_panels, aes(time, hr), colour = "black", linewidth = 1) +
  geom_line(data = res, aes(time, hr), colour = "red", linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = 1, colour = "grey60") +
  facet_wrap(~ method) +
  scale_x_continuous(breaks = seq(0, max(res$time), by = 24)) +
  labs(x = "Time (months)", y = "HR(Surgery)") +
  theme_bw()
save_figure(p, "Figure8_CRC_surgery_TVE.png")

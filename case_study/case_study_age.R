# CRC application, Model 1: non-linear effect (NLE) of age (Figure 7).
# Run from the repository root: Rscript case_study/case_study_age.R
suppressMessages({
  library(survival); library(rstpm2); library(imputer); library(mice); library(miceadds)
  library(ggplot2)
})
source("R/labels.R")
source("R/utils_case_study.R")

CRC <- read_CRC()
RNG_AGE <- range(CRC$Age)
DF_AGE <- 4    # df of the natural spline of age (see the AIC table in Appendix D)
AGE_REF <- 60  # reference age of the hazard ratio
AGE_GRID <- seq(max(25, RNG_AGE[1]), min(85, RNG_AGE[2]), by = 1)

#' Fit the NLE model and return HR(age vs AGE_REF) on AGE_GRID
#' (the HR does not depend on time in this model)
evaluate_age <- function(dat) {
  fit <- stpm2(Surv(time, event) ~ nsx(Age, df = DF_AGE, Boundary.knots = RNG_AGE) +
                 Surgery + MSI.MSS.status + Sex + Stage, data = dat)
  nd <- data.frame(time = 12, Age = AGE_GRID, Surgery = 0, MSI.MSS.status = 0,
                   Sex = 0, Stage = 1)
  ref <- transform(nd, Age = AGE_REF)
  data.frame(Age = AGE_GRID,
             hr = predict(fit, newdata = nd, type = "haz") /
               predict(fit, newdata = ref, type = "haz"))
}

set.seed(1234)
ampdat <- ampute_CRC(CRC, "Age")   # one amputed dataset, used by every method
reference <- evaluate_age(CRC)
res <- do.call(rbind, lapply(CASE_METHODS, impute_and_evaluate,
                             ampdat = ampdat, evaluate = evaluate_age))

# Mean absolute deviation from the reference curve on the log-HR scale
res$ref <- reference$hr[match(res$Age, reference$Age)]
mad <- aggregate(cbind(MAD = abs(log(hr) - log(ref))) ~ method, data = res, FUN = mean)
mad$method <- method_factor(mad$method)
mad <- mad[order(mad$MAD), ]
print(mad, digits = 2, row.names = FALSE)

write.csv(rbind(data.frame(reference, method = "complete_data"), res[, names(res) != "ref"]),
          results_path("tables", "case_study_age_curves.csv"), row.names = FALSE)
write.csv(mad, results_path("tables", "case_study_age_MAD.csv"), row.names = FALSE)

res$method <- method_factor(res$method)
ref_panels <- merge(reference, data.frame(method = levels(res$method)))
ref_panels$method <- factor(ref_panels$method, levels = levels(res$method))
p <- ggplot() +
  geom_line(data = ref_panels, aes(Age, hr), colour = "black", linewidth = 1) +
  geom_line(data = res, aes(Age, hr), colour = "red", linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = 1, colour = "grey60") +
  facet_wrap(~ method) +
  scale_y_log10() +
  labs(x = "Age (years)", y = sprintf("HR(Age vs %d)", AGE_REF)) +
  theme_bw()
save_figure(p, "Figure7_CRC_age_NLE.png")

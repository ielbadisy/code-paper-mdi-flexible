# Age effect in the CRC application fitted as a non-linear effect of age
# (natural spline of Age in the formula) and plotted against age.
# The AIC table compares linear Age with natural splines of Age (df 1 to 6).
suppressMessages({
  library(tidyverse); library(survival); library(rstpm2)
  library(imputer); library(mice); library(miceadds)
})
set.seed(1234)

CRC <- read.csv("CRC.csv") %>%
  imputer::imputer("complete") %>%
  select(time, event, Age, Surgery, MSI.MSS.status, Sex, Stage) %>%
  mutate(
    MSI.MSS.status = ifelse(MSI.MSS.status == "Yes", 1, 0),
    Stage = as.numeric(as.factor(Stage)),
    Surgery = ifelse(Surgery == "Yes", 1, 0),
    Sex = ifelse(Sex == "Male", 1, 0)
  )

# df for the spline of Age chosen by AIC on the full data
aic <- function(f) 2 * f@min + 2 * length(coef(f))
aic_tab <- map_dfr(1:6, function(d) {
  f <- stpm2(Surv(time, event) ~ nsx(Age, df = d, Boundary.knots = range(CRC$Age)) +
               Surgery + MSI.MSS.status + Sex + Stage, data = CRC)
  tibble(df_age = d, AIC = aic(f))
})
f_lin <- stpm2(Surv(time, event) ~ Age + Surgery + MSI.MSS.status + Sex + Stage,
               data = CRC)
aic_tab <- bind_rows(tibble(df_age = NA, AIC = aic(f_lin)), aic_tab)
aic_tab$model <- c("linear Age", sprintf("nsx(Age, df = %d)", 1:6))
print(as.data.frame(aic_tab[, c("model", "AIC")]), digits = 4)
DF_AGE <- 4   # df of the age spline

age_grid <- seq(25, 85, by = 1)
AGE_REF  <- 60
rng      <- range(CRC$Age)
age_grid <- age_grid[age_grid >= rng[1] & age_grid <= rng[2]]

# Model evaluation: HR(Age vs AGE_REF) across age (constant over time in this model)
evaluate_model <- function(impdat) {
  fit <- stpm2(Surv(time, event) ~ nsx(Age, df = DF_AGE, Boundary.knots = rng) +
                 Surgery + MSI.MSS.status + Sex + Stage, data = impdat)
  nd  <- data.frame(time = 12, Age = age_grid, Surgery = 0, MSI.MSS.status = 0, Sex = 0, Stage = 1)
  ref <- transform(nd, Age = AGE_REF)
  data.frame(Age = age_grid,
             hr_Age = predict(fit, newdata = nd, type = "haz") / predict(fit, newdata = ref, type = "haz"))
}

generate_NA_Age <- function(dat, mech = "MNAR", pmiss = 0.3, patterns = c(1, 1, 0, 1, 1, 1, 1)) {
  mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp
}

imputation_methods <- c("complete-cases", "naive", "hotdeck", "mice", "cart",
                        "famd", "knn", "misscforest", "missforest", "missranger",
                        "micecart", "micerf")
method_names <- c("complete-cases" = "CompleteCases", "naive" = "Mean/mode", "hotdeck" = "Hotdeck",
                  "mice" = "MICE", "cart" = "CART", "famd" = "FAMD", "knn" = "KNN",
                  "misscforest" = "missCforest", "missforest" = "missForest",
                  "missranger" = "missRanger", "micecart" = "miceCART", "micerf" = "miceRF")
facet_order <- unname(method_names)

process_imputation <- function(method) {
  if (method == "complete-cases") {
    impdat <- imputer::imputer(generate_NA_Age(CRC), "complete")
  } else if (method %in% c("micecart", "micerf", "mice")) {
    dat_mice <- mice::complete(imputer::imputer(generate_NA_Age(CRC), method), "all")
    result <- Reduce(`+`, lapply(dat_mice, evaluate_model)) / length(dat_mice)
    return(data.frame(result, method = method))
  } else {
    impdat <- imputer::imputer(generate_NA_Age(CRC), method)
  }
  data.frame(evaluate_model(impdat), method = method)
}

true_Age <- evaluate_model(CRC)
res_Age  <- map_dfr(imputation_methods, process_imputation) %>%
  mutate(method = factor(recode(method, !!!method_names), levels = facet_order))
true_exp <- true_Age %>% tidyr::crossing(method = facet_order) %>%
  mutate(method = factor(method, levels = facet_order))
write.csv(bind_rows(mutate(true_Age, method = "Reference"), res_Age),
          "res_Age_nsx.csv", row.names = FALSE)

ggplot() +
  geom_line(data = true_exp, aes(Age, hr_Age), colour = "black", linewidth = 1) +
  geom_line(data = res_Age, aes(Age, hr_Age), colour = "red", linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = 1, colour = "grey60") +
  facet_wrap(~ method) +
  scale_y_log10() +
  labs(x = "Age (years)", y = sprintf("HR(Age vs %d)", AGE_REF)) +
  theme_bw()
ggsave("fig7-crc-age-nle.png", width = 8, height = 6, dpi = 300)

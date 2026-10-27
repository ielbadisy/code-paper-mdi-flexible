# Load required packages
pacman::p_load(
  naniar, raster, png, gridExtra, ggplot2, flextable, tidyverse, finalfit, rmdHelpers, 
  texreg, tab, modelsummary, gtsummary, forestmodel, MatchIt, table1, survminer, 
  survival, broom, knitr, kableExtra, cowplot, glmpath, ggpubr, glmnet, missForest, 
  missRanger, randomForestSRC, ranger, rstpm2, imputer, miceadds, mice, purrr
)

# Set random seed for reproducibility
set.seed(1234)

# Load and preprocess data
CRC <- read.csv("CRC.csv") %>%
  imputer::imputer("complete") %>%
  select(time, event, Age, Surgery, MSI.MSS.status, Sex, Stage) %>%
  mutate(
    MSI.MSS.status = ifelse(MSI.MSS.status == "Yes", 1, 0),
    Stage = as.numeric(as.factor(Stage)),
    Surgery = ifelse(Surgery == "Yes", 1, 0),
    Sex = ifelse(Sex == "Male", 1, 0)
  )

# Define a function to evaluate models with TVE on Surgery
evaluate_tve_surgery <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, event) ~ Age + Surgery + MSI.MSS.status + Sex + Stage,
                       data = impdat, tvc = list(Surgery = 3))
  
  preds <- rstpm2::predict(fit, type = "hr",
                           newdata = data.frame(Age = 19, Surgery = 0, MSI.MSS.status = 0, Sex = 0, Stage = 1),
                           var = "Surgery", grid = TRUE, seqLength = 322, full = TRUE, se.fit = TRUE)
  
  data.frame(time_Surgery = preds$time, hr_Surgery = preds$Estimate)
}

# Generate missing data for Surgery
generate_NA_Surgery <- function(dat, mech = "MNAR", pmiss = 0.3, patterns = c(1, 1, 1, 0, 1, 1, 1)) {
  mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp
}

# List of imputation methods
imputation_methods <- c("complete-cases", "naive", "hotdeck", "mice", "cart",
                        "famd", "knn", "misscforest", "missforest", "missranger",
                        "micecart", "micerf")

# Function to process each imputation method
process_imputation_surgery <- function(method) {
  if (method == "complete-cases") {
    impdat <- imputer::imputer(generate_NA_Surgery(CRC), "complete")
  } else if (method %in% c("micecart", "micerf", "mice")) {
    mice_imp <- imputer::imputer(generate_NA_Surgery(CRC), method)
    dat_mice <- mice::complete(mice_imp, "all")
    result <- Reduce(`+`, lapply(dat_mice, evaluate_tve_surgery)) / length(dat_mice)
    return(data.frame(result, method = method))
  } else {
    impdat <- imputer::imputer(generate_NA_Surgery(CRC), method)
  }
  
  data.frame(evaluate_tve_surgery(impdat), method = method)
}

# Compute the TVE for the complete dataset
true_Surgery <- evaluate_tve_surgery(CRC)

# Run all imputation methods
res_Surgery <- map_dfr(imputation_methods, process_imputation_surgery)

# Standard method names for renaming
method_names <- c(
  "complete-cases" = "CompleteCases",
  "naive" = "Mean/mode",
  "hotdeck" = "Hotdeck",
  "mice" = "MICE",
  "cart" = "CART",
  "famd" = "FAMD",
  "knn" = "KNN",
  "misscforest" = "missCforest",
  "missforest" = "missForest",
  "missranger" = "missRanger",
  "micecart" = "miceCART",
  "micerf" = "miceRF"
)

# Define the correct facet order
facet_order <- c("CompleteCases", "Mean/mode", "Hotdeck", "MICE", "CART",
                 "FAMD", "KNN", "missCforest", "missForest", "missRanger",
                 "miceCART", "miceRF")

# Rename and order methods
res_Surgery_filtered <- res_Surgery %>%
  mutate(
    method = recode(method, !!!method_names),
    method = factor(method, levels = facet_order)
  )

true_Surgery_expanded <- true_Surgery %>%
  tidyr::crossing(method = facet_order) %>%
  mutate(method = factor(method, levels = facet_order))

# Merge datasets
res_Surgery_updated <- bind_rows(res_Surgery_filtered, true_Surgery_expanded)

# Final Plot
ggplot(res_Surgery_updated, aes(time_Surgery, hr_Surgery)) + 
  geom_smooth(data = true_Surgery_expanded, aes(time_Surgery, hr_Surgery),
              color = "black", linetype = "solid", size = 1, se = FALSE, span = 1) +
  geom_smooth(data = res_Surgery_filtered, aes(time_Surgery, hr_Surgery),
              color = "red", linetype = "dotted", se = FALSE, span = 1) +
  facet_wrap(~ method, scales = "free") +
  coord_cartesian(ylim = c(-0.5, 1.5), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, max(res_Surgery_updated$time_Surgery, na.rm = TRUE), by = 12)) +  # ✅ Add this line

  labs(x = "Time", y = "TVE(Surgery)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot
ggsave("final_NLE_Surgery_plot.png", width = 8, height = 6, dpi = 300)

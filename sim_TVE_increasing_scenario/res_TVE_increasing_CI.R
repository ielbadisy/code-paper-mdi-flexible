#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
TVE_data <- vroom("res_TVE_increasing.csv") %>%
  mutate(
    method = recode(method,
                    complete = "CompleteCases",
                    naive = "Mean/mode",
                    hotdeck = "Hotdeck",
                    mice = "MICE",
                    cart = "CART",
                    famd = "FAMD",
                    knn = "KNN",
                    misscforest = "missCforest",
                    missforest = "missForest",
                    missranger = "missRanger",
                    micecart = "miceCART",
                    micerf = "miceRF"),
    method = factor(method, levels = c("CompleteCases", "Mean/mode", "Hotdeck", "MICE", "CART",
                                       "FAMD", "KNN", "missCforest", "missForest", "missRanger",
                                       "miceCART", "miceRF"))
  ) %>%
  # apply the correct formula for increasing shape
  mutate(
    true_hr_x3 = exp(0.01 + 0.65 * time^0.4),    
    # compute coverage correctly for CompleteCases only
    coverage_hr_x3 = if_else(
      method == "CompleteCases",
      between(true_hr_x3, hr_x3_lower, hr_x3_upper),
      coverage_hr_x3
    )
  )

#====================[ DEFINE TIMEPOINTS FOR MATCHING ]====================
time_points <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
tolerance <- 0.01

#====================[ PLOTTING FUNCTION ]====================
plot_coverage_tve <- function(mechanism) {
  
  coverage_results <- TVE_data %>%
    filter(mech == mechanism, method != "complete_data") %>%
    filter(sapply(time, function(t) any(abs(t - time_points) < tolerance))) %>%
    mutate(
      matched_time = sapply(time, function(t) time_points[which.min(abs(t - time_points))])
    ) %>%
    group_by(method, matched_time) %>%
    summarise(
      coverage_rate = mean(coverage_hr_x3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Time = matched_time, Coverage = coverage_rate)
  
  facet_plot <- ggplot(coverage_results, aes(x = Time, y = Coverage)) +
    geom_point(size = 2, shape = 17, color = "red") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    facet_wrap(~method, ncol = 4, scales = "fixed") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(
      breaks = seq(1, 5, by = 1),               # Show only full-year labels
      minor_breaks = seq(0.5, 5, by = 0.5)      # Optional: show minor tick marks
    ) +
    labs(
      x = "Time (years)",
      y = expression("Coverage rate of HR"[X[3]])
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  ggsave(paste0("coverage_TVE_increasing_", mechanism, ".png"),
         facet_plot, width = 8, height = 6, dpi = 300)
  
  message("CI coverage plot saved for TVE Increasing - ", mechanism)
}

#====================[ RUN FOR ALL MECHANISMS ]====================
c("MNAR", "MAR", "MCAR") %>% walk(plot_coverage_tve)

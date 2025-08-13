#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
TVE_data <- vroom("res_TVE_bathub.csv") %>%
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
  # Compute vector of true HR(X3) based on formula applied at each time
  mutate(
    true_hr_x3 = exp(0.01 + 1.01 * exp(-time / 0.3) + 0.06 * time^1.1),
    
    # Replace coverage only for CompleteCases using true HR at the same time point
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
  
  # Select and match data to time points
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
  
  # Plot
  facet_plot <- ggplot(coverage_results, aes(x = Time, y = Coverage)) +
    geom_point(size = 2, shape = 17, color = "red") +
    #geom_line(aes(group = method), color = "blue", linewidth = 0.4) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    facet_wrap(~method, ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(1, 5, by = 1)) +  # << Added this line
    labs(
      x = "Time (years)",
      y = expression("Coverage rate of HR"[X[3]])
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  ggsave(paste0("coverage_TVE_bathub_", mechanism, ".png"),
         facet_plot,  width = 8, height = 6, dpi = 300)
  
  message("CI coverage plot saved for TVE - ", mechanism)
}

#====================[ RUN FOR ALL MECHANISMS ]====================
c("MNAR", "MAR", "MCAR") %>% walk(plot_coverage_tve)

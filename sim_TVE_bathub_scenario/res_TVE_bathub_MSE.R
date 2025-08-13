#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
TVE_3030 <- vroom("res_TVE_bathub.csv") %>%
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
  )

#====================[ DEFINE TIMEPOINTS FOR MATCHING ]====================
time_points <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
tolerance <- 0.01

#====================[ PLOTTING FUNCTION ]====================
plot_mse <- function(mechanism) {
  
  mse_data <- TVE_3030 %>%
    filter(mech == mechanism, !is.na(method)) %>%
    mutate(
      # theoretical HR for bathtub-shaped time-varying effect
      true_hr = exp(0.01 + 1.01 * exp(-time_x3 / 0.3) + 0.06 * (time_x3 ^ 1.1)),
      
      # use true HR for all methods including CompleteCases
      squared_error = (hr_x3 - true_hr)^2,
      
      # match each time point to the closest predefined grid value
      matched_time = sapply(time_x3, function(t) time_points[which.min(abs(t - time_points))])
    ) %>%
    filter(abs(time_x3 - matched_time) < tolerance)
  
  #====================[ BOXPLOT ]====================
  mse_boxplot <- ggplot(mse_data, aes(x = matched_time, y = squared_error, group = matched_time)) +
    geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.size = 0.1) +
    facet_wrap(~method, ncol = 4) +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(breaks = seq(1, 5, by = 1)) +  # Show only full years
    labs(
      x = "Time (years)",
      y = expression(paste("MSE of ", HR[X[3]]))
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  ggsave(paste0("MSE_TVE_bathub_", mechanism, "_boxplot.png"),
         mse_boxplot, width = 8, height = 6, dpi = 300)
  
  message("MSE grid-based plot saved for TVE Bathtub - ", mechanism)
}

#====================[ RUN FOR ALL MECHANISMS ]====================
c("MCAR", "MAR", "MNAR") %>% walk(plot_mse)
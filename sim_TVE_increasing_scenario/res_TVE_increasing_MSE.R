#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
TVE_3030 <- vroom("res_TVE_increasing.csv") %>%
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

#====================[ DEFINE TIMEPOINTS ]====================
time_points <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
tolerance <- 0.01

#====================[ PLOTTING FUNCTION ]====================
plot_mse <- function(mechanism) {
  
  mse_data <- TVE_3030 %>%
    filter(mech == mechanism, !is.na(method)) %>%
    filter(sapply(time_x3, function(t) any(abs(t - time_points) < tolerance))) %>%
    mutate(
      matched_time = sapply(time_x3, function(t) time_points[which.min(abs(t - time_points))]),
      
      # True HR only for CompleteCases using correct decreasing shape formula
      true_hr = exp(0.01 + 0.65 * matched_time^0.4),
      
      squared_error = if_else(
        method == "CompleteCases",
        (hr_x3 - true_hr)^2,
        (hr_x3 - obs_hr_x3)^2
      )
    )
  
  mse_boxplot <- ggplot(mse_data, aes(x = matched_time, y = squared_error, group = matched_time)) +
    geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.size = 0.2) +
    facet_wrap(~method, ncol = 4) +
    scale_y_continuous(limits = c(0, 3)) +
    scale_x_continuous(
      breaks = seq(1, 5, by = 1),               # ✅ Show integer ticks
      minor_breaks = seq(0.5, 5, by = 0.5)      # ✅ Optional minor ticks at 0.5
    ) +
    labs(
      x = "Time (years)",
      y = expression(paste("MSE of ", HR[X[3]]))
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  ggsave(paste0("MSE_TVE_increasing_", mechanism, "_boxplot.png"),
         mse_boxplot, width = 8, height = 6, dpi = 300)
  
  message("MSE plot saved for TVE  Increasing - ", mechanism)
}

#====================[ RUN FOR ALL MECHANISMS ]=================================
c("MCAR", "MAR", "MNAR") %>% walk(plot_mse)

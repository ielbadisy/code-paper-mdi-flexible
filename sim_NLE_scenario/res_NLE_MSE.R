#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
NLE_3030 <- vroom("res_NLE.csv") %>%
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

#====================[ DEFINE MATCHING GRID FOR X2 ]====================
x2_points <- -4:6
tolerance <- 0.01

#====================[ PLOTTING FUNCTION ]====================
plot_mse_NLE <- function(mechanism) {
  
  mse_data <- NLE_3030 %>%
    filter(mech == mechanism, method != "complete_data") %>%
    mutate(
      # Compute true HR for x2b (based on time_x2)
      true_hr_x2 = 0.01 + 1.01 * exp(-time_x2 / 0.3) + 0.06 * (time_x2 ^ 1.3),
      
      # Compute squared error
      squared_error = (hr_x2 - true_hr_x2)^2,
      
      # Match x2b to nearest point on the grid
      matched_x2 = sapply(x2b, function(x) x2_points[which.min(abs(x - x2_points))])
    ) %>%
    filter(abs(x2b - matched_x2) < tolerance)
  
  #====================[ BOXPLOT ]====================
  mse_boxplot <- ggplot(mse_data, aes(x = as.factor(matched_x2), y = squared_error)) +
    geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.size = 0.1) +
    facet_wrap(~method, ncol = 4) +
    scale_y_continuous(limits = c(0, 3)) +
    labs(
      x = expression(X[2]),
      y = expression(paste("MSE of ", HR[X[2]]))
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  # Save the plot
  ggsave(paste0("MSE_NLE_x2_grid_", mechanism, ".png"),
         mse_boxplot, width = 8, height = 6, dpi = 300)
  
  message("MSE grid-based plot saved for ", mechanism)
}

#====================[ RUN FOR ALL MECHANISMS ]====================
c("MCAR", "MAR", "MNAR") %>% walk(plot_mse_NLE)

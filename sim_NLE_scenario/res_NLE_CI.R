#====================[ LOAD REQUIRED PACKAGES ]====================
pacman::p_load(tidyverse, vroom)

#====================[ LOAD AND FORMAT DATASET ]====================
# Load main results
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
                                       "miceCART", "miceRF")),
    x2b_round = round(x2b, 2)
  )

# Load true HR for CompleteCases only
original_truth <- vroom("original_NLE.csv") %>%
  as_tibble() %>%
  mutate(x2b_round = round(x2b, 2)) %>%
  rename(true_hr_x2_completecases = hr_x2) %>%
  dplyr::select(x2b_round, true_hr_x2_completecases)

#====================[ DEFINE X2 BINNING STRATEGY ]====================
breaks_x2 <- c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6)
labels_x2 <- c("[-3,-2)", "[-2,-1)", "[-1,0)", "[0,1)", "[1,2)", 
               "[2,3)", "[3,4)", "[4,5)", "[5,6]")

#====================[ PLOTTING FUNCTION ]====================
plot_coverage_x2_discrete <- function(mechanism, centers = -3:6, tol = 0.1) {
  
  # Merge with truth and compute logical coverage
  NLE_adj <- NLE_3030 %>%
    filter(mech == mechanism) %>%
    left_join(original_truth, by = "x2b_round") %>%
    mutate(
      coverage_final = as.numeric(if_else(
        method == "CompleteCases",
        true_hr_x2_completecases >= hr_x2_lower & true_hr_x2_completecases <= hr_x2_upper,
        coverage_hr_x2
      ))
    ) %>%
    filter(!is.na(x2b), !is.na(method))
  
  # Compute coverage within ±tol around each integer center
  coverage_results <- map_dfr(centers, function(center) {
    NLE_adj %>%
      filter(x2b >= center - tol, x2b < center + tol) %>%
      group_by(method) %>%
      summarise(
        x2b_center = center,
        coverage_rate = mean(coverage_final, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  })
  
  # Plot
  facet_plot <- ggplot(coverage_results, aes(x = factor(x2b_center), y = coverage_rate)) +
    geom_point(size = 2, shape = 17, color = "red") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    facet_wrap(~method, ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(X[2]),
      y = expression("Coverage rate of HR"[X[2]])
    ) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  ggsave(paste0("coverage_NLE_x2_integer_", mechanism, ".png"),
         facet_plot, width = 8, height = 6, dpi = 300)
  
  message("[OK] CI coverage plot by integer X2 (±", tol, ") saved for ", mechanism)
}

#====================[ GENERATE PLOTS ]====================

c("MCAR", "MAR", "MNAR") %>% walk(~plot_coverage_x2_discrete(.x, tol = 0.1))

pacman::p_load(tidyverse, vroom)

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

#=========================== setting 1: PH ==========================
PH_3030 <- vroom("res_PH.csv")
PH_3030 <- PH_3030 %>% mutate(covariates = recode(covariates, 
                                      x1 = "X1",
                                      x2 = "X2 continuous",
                                      x31 = "X3 binary"))

PH_3030 %>%
  filter(covariates != "x1") -> PH_3030

PH_3030 <- PH_3030 %>% mutate(method = recode(method, 
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
                                      micerf = "miceRF"))

PH_3030$method <- factor(PH_3030$method, levels=c("CompleteCases", "Mean/mode", "Hotdeck", "MICE", "CART", 
                                          "FAMD", "KNN", "missCforest",  "missForest", "missRanger",
                                          "miceCART", "miceRF"))

## MCAR -----------------------------------------------------------------------------------------------
# Bias plot
Long_res_bias <- gather(PH_3030, metric, value, c(bias))
ggplot(Long_res_bias, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  geom_vline(lty = 20, aes(xintercept = 0), colour = 'red') +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_Bias.png", width = 8, height = 6, dpi = 300)

# Empirical Error plot
Long_res_se_est <- gather(PH_3030, metric, value, c(se_est))
ggplot(Long_res_se_est, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  coord_cartesian(xlim = c(0, 0.15), expand = FALSE) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_Empirical_Error.png", width = 8, height = 6, dpi = 300)

# MSE plot
Long_res_MSE_estimate <- gather(PH_3030, metric, value, c(MSE_estimate))
ggplot(Long_res_MSE_estimate, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  geom_vline(lty = 20, aes(xintercept = 0), colour = 'red') +
  coord_cartesian(xlim = c(0, 0.3), expand = FALSE) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_MSE.png", width = 8, height = 6, dpi = 300)

# IBS plot
Long_res_IBS <- gather(PH_3030, metric, value, c(IBS))
ggplot(Long_res_IBS, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_IBS.pNG", width = 8, height = 6, dpi = 300)

# AUC plot
Long_res_AUC <- gather(PH_3030, metric, value, c(AUC))
ggplot(Long_res_AUC, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_AUC.png", width = 8, height = 6, dpi = 300)

# C-index plot
Long_res_cindex <- gather(PH_3030, metric, value, c(cindex))
ggplot(Long_res_cindex, aes(x = value, group = factor(method))) +
  geom_boxplot(aes(fill = method), outlier.shape = NA) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_grid(mech ~ covariates) +
  theme_bw()+
  theme(axis.text.y = element_blank()) 
ggsave("PH_Cindex.png", width = 8, height = 6, dpi = 300)


### All the results in a table

#### Table 1
PH_3030 |> group_by(covariates, method, mech) |>
filter(covariates == "x2_continuous") |>
   summarise(Bias = mean(bias),
             Empirical_SE = mean(se_est),
             Relative_bias = mean(bias_rel),
             MSE = mean(MSE_estimate),
             CI_length = mean(lo95 - hi95),
             CI_coverage = mean(cover),
             IBS = mean(IBS),
             AUC = mean(AUC),
             Cindex = mean(cindex)
             ) |>
  round_df(4)

#### Table 2
PH_3030 |> group_by(covariates, method, mech) |>
  filter(covariates == "x3_binary") |>
  summarise(Bias = mean(bias),
            Empirical_SE = mean(se_est),
            Relative_bias = mean(bias_rel),
            MSE = mean(MSE_estimate),
            CI_length = mean(lo95 - hi95),
            CI_coverage = mean(cover),
            IBS = mean(IBS),
            AUC = mean(AUC),
            Cindex = mean(cindex)
            ) |>
  round_df(4)
  



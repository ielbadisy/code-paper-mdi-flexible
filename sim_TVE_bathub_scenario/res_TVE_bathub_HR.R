pacman::p_load(tidyverse, vroom)
#=========================== setting 3: TVE ==========================
TVE_3030 <- vroom("res_TVE_bathub.csv")

TVE_3030 <- TVE_3030 %>% mutate(method = recode(method,
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

TVE_3030$method <- factor(TVE_3030$method, levels=c("CompleteCases", "Mean/mode", "Hotdeck", "MICE", "CART", 
                                                    "FAMD", "KNN", "missCforest",  "missForest", "missRanger",
                                                    "miceCART", "miceRF"))

### MCAR -------------------------------------------------------------------
TVE_3030 %>%
  subset(mech == "MCAR" & method != "complete_data") %>%
  mutate(avg = approx(time_x3, hr_x3,xout=time_x3)$y) %>%
  ggplot(aes(time_x3, hr_x3))+ 
  geom_line(aes(time_x3, hr_x3, group=run_id), colour = "grey")+
  geom_smooth(se=FALSE, colour="red", linetype="11")+
  facet_wrap(facets = vars(method))+
  stat_function(fun = function(x) {exp(0.01 + 1.01 * exp(-x/0.3) + 0.06 * (x^1.1))}) +
  coord_cartesian(ylim=c(0, 5), expand = TRUE)+
  labs(x = "Time (years)",
       y = expression(HR[X[3]]))+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
ggsave("TVE_MCAR_bathub.png", width = 8, height = 6, dpi = 300)

### MAR -------------------------------------------------------------------
TVE_3030 %>%
  subset(mech == "MAR" & method != "complete_data") %>%
  mutate(avg = approx(time_x3, hr_x3,xout=time_x3)$y) %>%
  ggplot(aes(time_x3, hr_x3))+ 
  geom_line(aes(time_x3, hr_x3, group=run_id), colour = "grey")+
  geom_smooth(se=FALSE, colour="red", linetype="11")+
  facet_wrap(facets = vars(method))+
  stat_function(fun = function(x) {exp(0.01 + 1.01 * exp(-x/0.3) + 0.06 * (x^1.1))}) +
  coord_cartesian(ylim=c(0, 5), expand = TRUE)+
  labs(x = "Time (years)",
       y = expression(HR[X[3]]))+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
ggsave("TVE_MAR_bathub.png", width = 8, height = 6, dpi = 300)

### MNAR -------------------------------------------------------------------
TVE_3030 %>%
  subset(mech == "MNAR" & method != "complete_data") %>%
  mutate(avg = approx(time_x3, hr_x3,xout=time_x3)$y) %>%
  ggplot(aes(time_x3, hr_x3))+ 
  geom_line(aes(time_x3, hr_x3, group=run_id), colour = "grey")+
  geom_smooth(se=FALSE, colour="red", linetype="11")+
  facet_wrap(facets = vars(method))+
  stat_function(fun = function(x) {exp(0.01 + 1.01 * exp(-x/0.3) + 0.06 * (x^1.1))}) +
  coord_cartesian(ylim=c(0, 5), expand = TRUE)+
  labs(x = "Time (years)",
       y = expression(HR[X[3]]))+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
ggsave("TVE_MNAR_bathub.png", width = 8, height = 6, dpi = 300)
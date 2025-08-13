pacman::p_load(tidyverse, vroom)
#=========================== setting 3: NLE ==========================
NLE_3030 <- vroom("res_NLE.csv")
original_NLE <- read.csv("original_NLE.csv")


NLE_3030 <- NLE_3030 %>% mutate(method = recode(method,
                                                complete_data = "CompleteData",
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

NLE_3030$method <- factor(NLE_3030$method, levels=c("CompleteData", "CompleteCases", "Mean/mode", "Hotdeck", "MICE", "CART", 
                                                    "FAMD", "KNN", "missCforest",  "missForest", "missRanger",
                                                    "miceCART", "miceRF"))

### MCAR -------------------------------------------------------------------
NLE_3030 %>%
  subset(mech == "MCAR" & method != "CompleteData") %>%
  ggplot(aes(x2b, hr_x2))+ 
   geom_line(aes(x2b, hr_x2, group=run_id), colour = "grey")+ # , color = method : changer the color
   geom_line(aes(x2b, hr_x2), original_NLE)+
   facet_wrap(facets = vars(method))+
   geom_smooth(se=FALSE, colour="red", span=0.3, linetype="11")+
   coord_cartesian(ylim=c(1, 3.5), expand = TRUE)+
   labs(
     x =expression(X[2]), 
     y = expression(HR[X[2]])
   )+ 
   theme(plot.title = element_text(hjust = 0.5))+
   theme_bw()
ggsave("NLE_MCAR.png", width = 8, height = 6, dpi = 300)

### MAR -------
NLE_3030 %>%
  subset(mech == "MAR" & method != "CompleteData") %>%
  ggplot(aes(x2b, hr_x2))+ 
  geom_line(aes(x2b, hr_x2, group=run_id), colour = "grey")+ # , color = method : changer the color
  geom_line(aes(x2b, hr_x2), original_NLE)+
  facet_wrap(facets = vars(method))+
  geom_smooth(se=FALSE, colour="red", span=0.3, linetype="11")+
  coord_cartesian(ylim=c(1, 3.5), expand = TRUE)+
  labs(
    x =expression(X[2]), 
    y = expression(HR[X[2]])
  )+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
ggsave("NLE_MAR.png", width = 8, height = 6, dpi = 300)


### MNAR ------

NLE_3030 %>%
  subset(mech == "MNAR" & method != "CompleteData") %>%
  ggplot(aes(x2b, hr_x2))+ 
  geom_line(aes(x2b, hr_x2, group=run_id), colour = "grey")+ # , color = method : changer the color
  geom_line(aes(x2b, hr_x2), original_NLE)+
  facet_wrap(facets = vars(method))+
  geom_smooth(se=FALSE, colour="red", span=0.3, linetype="11")+
  coord_cartesian(ylim=c(1, 3.5), expand = TRUE)+
  labs(
    x =expression(X[2]), 
    y = expression(HR[X[2]])
  )+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()
ggsave("NLE_MNAR.png", width = 8, height = 6, dpi = 300)

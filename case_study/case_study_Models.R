
## Models 
pacman::p_load(tidyverse)
set.seed(1234)
CRC <- read.csv("CRC.csv")
CRC <- imputer::imputer(CRC, "complete")
CRC <- select(CRC, time, event, Age, Surgery, MSI.MSS.status, Sex, Stage) 
CRC$MSI.MSS.status <- ifelse(CRC$MSI.MSS.status == "Yes", 1, 0)
CRC$Stage <- as.numeric(as.factor(CRC$Stage))
CRC$Surgery <- ifelse(CRC$Surgery=="Yes",1,0)
CRC$Sex <- ifelse(CRC$Sex=="Male",1,0)



### Model 1: NLE Age 

fit_Age <- rstpm2::stpm2(Surv(time, event) ~  Age + Surgery + MSI.MSS.status + Sex + Stage,
                         data = CRC,
                         tvc = list(Age = 4))

broom::tidy(fit_Age)




### Model 2: TVE 

fit_Surgery <- rstpm2::stpm2(Surv(time, event) ~  Age + Surgery + MSI.MSS.status + Sex + Stage,
                             data = CRC,
                             tvc = list(Surgery = 3))
broom::tidy(fit_Surgery)



# utils last edit: 07/02/2024
#---------- generate complete data -------------
genPH <- function(n = 1000, maxt = 5, hr = c(0.8, 0.5, 0.2), cens = 0.3){
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5) # chemotherapy

  delta=0.001
  t.grid=seq(0,maxt-delta,delta) #fine grid of times used to generate event times in a piecewise fashion

  # control censoring input values
  stopifnot("Only three levels of censoring are supported : 0.1 = 10%, 0.3 = 30% and 0.5 = 50%" = cens %in% c(0.1, 0.3, 0.5))

  # adjusments for obtaining the required censoring rate | lambda.c = rate of exponential censoring dist | lambda = baseline
  if(cens == 0.1){lambda.c <- 0.001 ; lambda = 0.4} # ~ 10%
  if(cens == 0.3){lambda.c <- 0.08 ; lambda = 0.2} # ~ 30%
  if(cens == 0.5){lambda.c <- 0.2 ; lambda = 0.1} # ~ 50%

  suppressWarnings({
    t.event.pw=matrix(nrow=n,ncol=length(t.grid))
    for(i in 1:length(t.grid)){
      u=runif(n,0,1)
      t.event.pw.a <- t.grid[i]-(log(u)/(lambda*exp(hr[1]*x1 + hr[2]*x2 + hr[3]*x3)))
      t.event.pw[,i]=ifelse(t.event.pw.a>=t.grid[i] & t.event.pw.a<(t.grid[i]+delta),t.event.pw.a,NA)
    }
    t.event=apply(t.event.pw,1,function(x) min(x,na.rm=T))
    t.event=ifelse(t.event==Inf,10,t.event)

    u.cens=runif(n,0,1)
    t.cens=-log(u.cens)/lambda.c

    d.a=ifelse(t.event<t.cens,1,2)
    d.a=ifelse(pmin(t.event,t.cens)>maxt,maxt,d.a)
    d=ifelse(d.a==1,1,0)

    t=pmin(t.event,t.cens)
    t=ifelse(t>maxt,maxt,t)

        dat=data.frame(time = t,status =d,x1,x2, x3)
    dat$cumhaz=mice::nelsonaalen(dat, time, status)
  })
  dat
}
#---------- introduce NA ------------
generate_NA_PH <- function(dat, mech, pmiss, patterns = c(1, 1, 1, 0, 0, 1)) {
  suppressWarnings({
    m = c("MAR", "MCAR", "MNAR")
    stopifnot(mech %in% m)
    if (mech == "MAR"){
      dat <- mice::ampute(dat, mech = "MAR", prop = pmiss, patterns = patterns)$amp
      ampdat <- dat
      return(ampdat)
    }
    if (mech == "MCAR"){
      dat <- mice::ampute(dat, mech = "MCAR", prop = pmiss, patterns = patterns)$amp
      ampdat <- dat
      return(ampdat)
    }
    if (mech == "MNAR"){
      dat <- mice::ampute(dat, mech = "MNAR", prop = pmiss, patterns = patterns)$amp
      ampdat <- dat
      return(ampdat)
    }
  })
}

#------------- estimation
## PH
evaluate_coxest_PH <- function(data, truelogHR){
  myformula <- as.formula(survival::Surv(time, status) ~ x1 + x2 + x3)
  dat <- data
  #need to check for arguments!!

  # Full data analysis
  coefs <- as.data.frame(summary(survival::coxph(myformula, data = dat))$coef)

  # return a data.frame of coefficients (est), upper and lower 95% limits
  out <- data.frame(covariates = row.names(coefs),
                    est = coefs$coef,
                    se_est = coefs$`se(coef)`,
                    lo95 = (coefs$coef + qnorm(0.025) * coefs$`se(coef)`),
                    hi95 = (coefs$coef + qnorm(0.975) * coefs$`se(coef)`)
  )
  out$bias = coefs$coef - truelogHR
  out$bias_rel =  ((coefs$coef - truelogHR) / truelogHR) * 100
  out$MSE_estimate = (coefs$coef - truelogHR)^2
  out$cover <- truelogHR >= out$lo95 & truelogHR <= out$hi95
  data.frame(out)
}
evaluate_coxest_PH_multiple <- function(mids, truelogHR){
  # pooled results
  fit <- with(mids, coxph(Surv(time, status) ~ x1 + x2 + x3))
  res <- summary(pool(fit))

  # return a data.frame of coefficients (est), upper and lower 95% limits
  out <- data.frame(covariates = as.character(res$term),
                    est = res$estimate,
                    se_est = res$std.error,
                    lo95 = (res$estimate + qnorm(0.025) * res$std.error),
                    hi95 = (res$estimate + qnorm(0.975) * res$std.error)
  )
  out$bias = out$est - truelogHR
  out$bias_rel =  ((out$est - truelogHR) / truelogHR)*100
  out$MSE_estimate = (out$est - truelogHR)^2
  out$cover <- truelogHR >= out$lo95 & truelogHR <= out$hi95
  data.frame(out)
}

#--------------- prediction
### PH
evaluate_coxperf_PH <- function(dat){
  myformula <- as.formula(survival::Surv(time, status) ~ x1 + x2 + x3)
  N = nrow(dat)

  #Initialization
  index.train = sample(1:N,2/3*N)
  data.train = dat[index.train,]
  data.test = dat[-index.train,]
  dis_time = sort(data.train$time[data.train$status == 1])  #the default time points

  ##---IBS
  # Gerds, T. A. and M. Schumacher (2006). Consistent estimation of the expected Brier score in general survival models with right-censored event times.
  # Biometrical Journal 48, 1029–1040.
  fitcox <- survival::coxph(myformula, data = data.train, x = TRUE)
  predcox <- predict(fitcox)
  predcoxnew <- predict(fitcox, newdata=data.test)
  surv_obj <- survival::Surv(data.train$time, data.train$status)
  surv_obj_new <- survival::Surv(data.test$time, data.test$status)
  IBS = survAUC::predErr(surv_obj, surv_obj_new, predcox, predcoxnew, dis_time, type = "brier", int.type = "weighted")$ierror

  ##---survAUC::AUC.cd()
  ### Chambless, L. E. and G. Diao (2006). Estimation of time-dependent area under the ROC curve for long-term risk prediction. Statistics in Medicine 25, 3474–3486
  fitcox = survival::coxph(myformula, data = data.train, x = TRUE)
  predcox <- predict(fitcox)
  predcoxnew <- predict(fitcox, newdata=data.test)
  surv_obj <- survival::Surv(data.train$time, data.train$status)
  surv_obj_new <- survival::Surv(data.test$time, data.test$status)
  AUC = survAUC::AUC.cd(surv_obj, surv_obj_new, predcox, predcoxnew, dis_time)$iauc

  ##---survAUC::preErr (absolute deviation between predicted and observed survival)
  # Schmid, M., T. Hielscher, T. Augustin, and O. Gefeller (2011).
  # implemented by hand see example here : https://stackoverflow.com/questions/44738753/obtaining-absolute-deviation-from-mean-for-two-sets-of-scores
  # A robust alter- native to the Schemper-Henderson estimator of prediction error. Biometrics 67, 524–535.
  fitcox <- survival::coxph(myformula, data = data.train, x = TRUE)
  predcox <- predict(fitcox)
  predcoxnew <- predict(fitcox, newdata=data.test)
  surv_obj <- survival::Surv(data.train$time, data.train$status)
  surv_obj_new <- survival::Surv(data.test$time, data.test$status)
  predError = survAUC::predErr(surv_obj, surv_obj_new, predcox, predcoxnew, dis_time, type = "robust", int.type = "unweighted")$ierror
  data.frame(IBS,
             AUC,
             predError)
}
extract_cindex_PH <- function(dat){
  cindex <- (survival::coxph(Surv(time, status) ~ x1 + x2 + x3, dat))$concordance[6]
  data.frame(cindex)
}

# utils last edit: 07/02/2024
#---------- generate complete data -------------
genTVE <- function(n = 1000, maxt = 5, shape = shape, cens = 0.3){
  x1 <- rnorm(n)
  x2 <- x1^2 + x1 + runif(n)
  x3 <- rbinom(n, 1, 0.5) # chemotherapy
  
  delta=0.01
  t.grid=seq(0,maxt-delta,delta) #fine grid of times used to generate event times in a piecewise fashion
  
  # control censoring input values
  stopifnot("Only three levels of censoring are supported : 0.1 = 10%, 0.3 = 30% and 0.5 = 50%" = cens %in% c(0.1, 0.3, 0.5))
  # adjusments for obtaining the required censoring rate | lambda.c = rate of exponential censoring dist | lambda = baseline
  if(cens == 0.1){lambda.c <- 0.02 ; lambda = 0.4} # ~ 10%
  if(cens == 0.3){lambda.c <- 0.15 ; lambda = 0.2} # ~ 30%
  if(cens == 0.5){lambda.c <- 0.2 ; lambda = 0.1} # ~ 50%
  
  suppressWarnings({
    t.event.pw=matrix(nrow=n,ncol=length(t.grid))
    for(i in 1:length(t.grid)){
      u=runif(n,0,1)
      t.event.pw.a <- switch(as.character(shape),
                             bathub = t.grid[i]-(log(u)/(lambda*exp(0.8*x1 + 0.5*x2 + (0.01+1.01*exp(-t.grid[i]/0.3)+0.06*(t.grid[i]^1.1))*x3))),

                             # alternative : curve(exp((0.01+1.01*exp(-x/0.1)+0.06*(x^1.1))), 0.1, to = 5, add = TRUE, ylim = c(0, 4), lty = 3, lwd = 2, col = "blue")
                             #bathtub: curve(exp((0.01+1.01*exp(-x/0.3)+0.06*(x^1.5))), 0.1, to = 5, add = TRUE, ylim = c(0, 4), lty = 3, lwd = 2, col = "blue")
                             decreasing = t.grid[i]-log(u)/(lambda*exp(0.8*x1 + 0.5*x2 + 0.01+(1.2/exp(t.grid[i])-0.02*t.grid[i]^0.7)*x3)),
                             #decreasing :  curve(exp(0.01+(1.2/exp(x-0.02*x^0.7))), 0.1, to = 5, add = TRUE, ylim = c(0, 4), lty = 3, lwd = 2, col = "blue")
                             increasing = t.grid[i]-log(u)/(lambda*exp(0.8*x1 + 0.5*x2 + (0.01+0.65*t.grid[i]^0.4)*x3))
                             # increasing : curve(exp(0.01+(0.65*x^0.4)), 0.1, to = 5, add = TRUE, ylim = c(0, 4), lty = 3, lwd = 2, col = "blue")
      )
      t.event.pw[,i]=ifelse(t.event.pw.a>=t.grid[i] & t.event.pw.a<(t.grid[i]+delta),t.event.pw.a,NA)
    }
    
    t.event=apply(t.event.pw,1,function(x) min(x,na.rm=T))
    t.event=ifelse(t.event==Inf, 30,t.event) # 10 could be a hypothetical  maximum observation time, its value have a direct impact on the b0
    
    u.cens=runif(n,0,1)
    t.cens=-log(u.cens)/lambda.c
    
    d.a=ifelse(t.event<t.cens,1,2)
    d.a=ifelse(pmin(t.event,t.cens)>maxt,maxt,d.a)
    d=ifelse(d.a==1,1,0)
    
    t=pmin(t.event,t.cens)
    t=ifelse(t>maxt,maxt,t)
    
    dat=data.frame(time = t,status =d, x1,x2, x3)
    cumhaz=mice::nelsonaalen(dat, time, status)
    cumhaz=ifelse(is.na(cumhaz), 0,cumhaz)
    dat$cumhaz <- cumhaz
  })
  dat
}


#---------- introduce NA ------------
generate_NA_TVE <- function(dat, mech, pmiss, patterns = c(1, 1, 1, 1, 0, 1)) {
  suppressWarnings({
    m = c("MAR", "MCAR", "MNAR")
    stopifnot(mech %in% m)
    # Set the missing data mechanism dynamically
    ampdat <- mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp
    return(ampdat)
  })
}

#------------- estimation
## TVE
evaluate_coxest_TVE <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, status) ~ x1 + x2 + x3,
                       data = impdat,
                       tvc = list(x3 = 7))
  h0 <- rstpm2::predict(fit, type="haz", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0),  grid=TRUE,
                        full=TRUE, se.fit=TRUE)$Estimate
  time <- rstpm2::predict(fit, type="haz", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0),  grid=TRUE,
                          full=TRUE, se.fit=TRUE)$time
  
  
  hr_x3 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                           full=TRUE, se.fit=TRUE)$Estimate
  hr_x3_lower <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                                 full=TRUE, se.fit=TRUE)$lower
  
  hr_x3_upper <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                                 full=TRUE, se.fit=TRUE)$upper
  time_x3 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                             full=TRUE, se.fit=TRUE)$time
  out <- cbind(time, h0, hr_x3, hr_x3_lower, hr_x3_upper, time_x3)
  
  return(out)
}

#------------ estimation fully observed covariate 
obs_hr_x3 <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, status) ~ x1 + x2 + x3,
                       data = impdat,
                       tvc = list(x3 = 7))
  
  
  obs_hr_x3 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                               full=TRUE, se.fit=TRUE)$Estimate
  obs_hr_x3_lower <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                                     full=TRUE, se.fit=TRUE)$lower
  
  obs_hr_x3_upper <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x3",  grid=TRUE,
                                     full=TRUE, se.fit=TRUE)$upper
  out <- cbind(obs_hr_x3, obs_hr_x3_lower, obs_hr_x3_upper)
  
  return(out)
}

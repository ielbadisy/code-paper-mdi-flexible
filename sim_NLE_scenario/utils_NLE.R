#---------- generate complete data -------------
genNLE <- function(n = 600, maxt = 5, cens = 0.3){
  
  x1 <- rnorm(n)
  x2 <- runif(n, -3, 6)
  x2 <- sort(x2)
  x3 <- rbinom(n, 1, 0.5) # chemotherapy
  
  delta=0.001
  t.grid=seq(0,maxt-delta,delta) #fine grid of times used to generate event times in a piecewise fashion
  
  # control censoring input values
  stopifnot("Only three levels of censoring are supported : 0.1 = 10%, 0.3 = 30% and 0.5 = 50%" = cens %in% c(0.1, 0.3, 0.5))
  
  # adjusments for obtaining the required censoring rate
  if(cens == 0.1){lambda.c <- 0.02 ; lambda = 0.4} # ~ 10%
  if(cens == 0.3){lambda.c <- 0.15 ; lambda = 0.2} # ~ 30%
  if(cens == 0.5){lambda.c <- 0.2 ; lambda = 0.1} # ~ 50%
  
  t.event.pw=matrix(nrow=n,ncol=length(t.grid))
  suppressWarnings({
    for(i in 1:length(t.grid)){
      u=runif(n,0,1)
      #t.event.pw.a <- t.grid[i]-(log(u)/(lambda*exp(0.8*x1 +(0.05+1.1*exp(-t.grid[i]/0.3)+0.04*(t.grid[i]^2))*x2 + 0.2*x3)))
      t.event.pw.a <- t.grid[i]-(log(u)/(lambda*exp(0.8*x1 + (0.01+1.01*exp(-t.grid[i]/0.3)+0.06*(t.grid[i]^1.3))*x2 + 0.3*x3)))
      t.event.pw[,i] <- ifelse(t.event.pw.a>=t.grid[i] & t.event.pw.a<(t.grid[i]+delta),t.event.pw.a,NA)
    }
    t.event=apply(t.event.pw,1,function(x) min(x,na.rm=T))
    t.event=ifelse(t.event==Inf,100,t.event)
    
    u.cens=runif(n,0,1)
    t.cens=-log(u.cens)/lambda.c
    
    d.a=ifelse(t.event<t.cens,1,2)
    d.a=ifelse(pmin(t.event,t.cens)>maxt,3,d.a)
    d=ifelse(d.a==1,1,0)
    
    t=pmin(t.event,t.cens)
    t=ifelse(t>maxt,maxt,t)
    
    dat=data.frame(time = t,status =d,x1,x2,x3)
    dat$cumhaz=mice::nelsonaalen(dat, time, status)
  })
  dat
}

#---------- introduce NA ------------

generate_NA_NLE <- function(dat, mech, pmiss, patterns = c(1, 1, 1, 0, 1, 1)) {
  suppressWarnings({
    m = c("MCAR", "MAR", "MNAR")
    stopifnot(mech %in% m)
    # Set the missing data mechanism dynamically
    ampdat <- mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp
    return(ampdat)
  })
}

#------------- estimation
evaluate_coxest_NLE <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, status) ~ x1 + x2 + x3,
                       data = impdat,
                       tvc = list(x2 = 7))
  
  n <- nrow(impdat) + 1
  
  h0 <- rstpm2::predict(fit, type="haz", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0),  grid=TRUE,seqLength=n,
                        full=TRUE, se.fit=TRUE)$Estimate
  time <- rstpm2::predict(fit, type="haz", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0),  grid=TRUE,seqLength=n,
                          full=TRUE, se.fit=TRUE)$time
  
  hr_x2 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                           full=TRUE, se.fit=TRUE)$Estimate
  hr_x2_lower <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                                 full=TRUE, se.fit=TRUE)$lower
  
  hr_x2_upper <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                                 full=TRUE, se.fit=TRUE)$upper
  
  time_x2 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                             full=TRUE, se.fit=TRUE)$time
  out <- cbind(time, h0, hr_x2, hr_x2_lower, hr_x2_upper, time_x2, x2b = sort(impdat$x2))
  return(out)
}

#------------ estimation fully observed covariate 

obs_hr_x2 <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, status) ~ x1 + x2 + x3,
                       data = impdat,
                       tvc = list(x2 = 7))
  n <- nrow(impdat) + 1
  
  obs_hr_x2 <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                               full=TRUE, se.fit=TRUE)$Estimate
  obs_hr_x2_lower <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                                     full=TRUE, se.fit=TRUE)$lower
  
  obs_hr_x2_upper <- rstpm2::predict(fit, type="hr", newdata=data.frame(x1 = 0, x2 = 0, x3 = 0), var = "x2",  grid=TRUE,seqLength=n,
                                     full=TRUE, se.fit=TRUE)$upper
  
  out <- cbind(obs_hr_x2, obs_hr_x2_lower, obs_hr_x2_upper)
  return(out)
}

#================ retreive coefs and se 

extract_coefs_NLE <- function(impdat) {
  fit <- rstpm2::stpm2(Surv(time, status) ~ x1 + x2 + x3,
                       data = impdat,
                       tvc = list(x2 = 7))
  beta=as.numeric(coef(fit))
  
  se=coef(summary(fit))
  
  out <- cbind(beta, se)
  return(out)
}

rm(list = ls())
source("utils_PH.R")
# import the utils
R = 2
N = 700
nCores = 32
miss = 0.3
cens = 0.3
mech = c("MCAR", "MAR", "MNAR")
set.seed = 15102024


pacman::p_load(dplyr, ggplot2, miceadds, rstpm2, stringr, biostat3, imputer, parallel, gower, fastDummies, mice, survival, mixgb, survAUC, purrr, pbapply, tictoc, missCforest)

#------------PH---------------
## single imputation ---------
grid_single_PH <- expand.grid(N = N,
                              mech = mech,
                              miss = miss,
                              cens = cens,
                              method = c("naive", "knn", "cart", "hotdeck", "famd", "missranger", "missforest", "misscforest"))

simcox_PH_single <- function(N, mech, miss, method, cens, run_id){
  truelogHR <- c(0.8, 0.5, 0.2)
  dat <- genPH(n = N, maxt = 5, cens = cens)
  ampdat <- generate_NA_PH(dat, mech = mech, pmiss = miss)
  ampdat$x3 <- as.factor(ampdat$x3)
  suppressWarnings({impdat <- imputer::imputer(ampdat, method = method)})
  res_est <- evaluate_coxest_PH(impdat, truelogHR)
  res_pred <- evaluate_coxperf_PH(impdat)
  cindex <- extract_cindex_PH(impdat)
  res_PH <- cbind(N, res_est, res_pred, cindex, method, miss, mech, cens, setting = "PH", run_id = run_id, row.names = NULL)
  return(res_PH)
}

grid_list <- pbapply::pblapply(seq_len(nrow(grid_single_PH)), function(i) {
  as.list(grid_single_PH[i, ])
})

simcox_PH_single2 <- purrr::possibly(.f = simcox_PH_single, otherwise = NULL)

out_PH_single <- pbapply::pblapply(1:R, function(run_id) {
  purrr::map_df(grid_list, ~ simcox_PH_single2(
    N = .x$N, 
    mech = .x$mech, 
    miss = .x$miss, 
    method = .x$method, 
    cens = .x$cens, 
    run_id = run_id  # pass the loop variable as run_id
  ))
}, cl = nCores)
# Combining the results into a single data frame
res_PH_single <- do.call(rbind, out_PH_single)

################################################################################

## complete cases imputation----------
grid_complete_PH <- expand.grid(N = N,
                                mech = mech,
                                miss = miss,
                                cens = cens,
                                method = c("complete"))

simcox_PH_complete <- function(N, mech, miss, method, cens, run_id){
  truelogHR <- c(0.8, 0.5, 0.2)
  dat <- genPH(n = N, maxt = 5, cens = cens)

  ampdat <- generate_NA_PH(dat, mech = mech, pmiss = miss)
  ampdat$x3 <- as.factor(ampdat$x3)
  
  suppressWarnings({impdat <- imputer::imputer(ampdat, method = method)})
  res_est <- evaluate_coxest_PH(impdat, truelogHR)
  res_pred <- evaluate_coxperf_PH(impdat)
  cindex <- extract_cindex_PH(impdat)
  res_PH <- cbind(N, res_est, res_pred, cindex, method, miss, mech, cens, setting = "PH", run_id = run_id, row.names = NULL)
  return(res_PH)
}

## simulation routine
simcox_PH_complete2 <- purrr::possibly(.f = simcox_PH_complete, otherwise = NULL)
grid_list <- pbapply::pblapply(seq_len(nrow(grid_complete_PH)), function(i) {
  as.list(grid_complete_PH[i, ])
})

# Now run the simulation using this list
out_PH_complete <- pbapply::pblapply(1:R, function(run_id) {
  purrr::map_df(grid_list, ~ simcox_PH_complete2(
    N = .x$N, 
    mech = .x$mech, 
    miss = .x$miss, 
    method = .x$method, 
    cens = .x$cens, 
    run_id = run_id  # pass the loop variable as run_id
  ))
})

res_PH_complete <- do.call(rbind, out_PH_complete)


################################################################################

## multiple imputation-----------
grid_multiple_PH <- expand.grid(N = N,
                                mech = mech,
                                miss = miss,
                                cens = cens,
                                method =  c("mice", "micecart", "micerf"))

simcox_PH_multiple <- function(N, mech, miss, method, cens, run_id){
  truelogHR <- c(0.8, 0.5, 0.2)
  dat <- genPH(n = N, maxt = 5, cens = cens)
  ampdat <- generate_NA_PH(dat, mech = mech, pmiss = miss)
  ampdat$x3 <- as.factor(ampdat$x3)
  
  suppressWarnings({impdat_mids <- imputer::imputer(ampdat, method = method)})
  impdat_list <- miceadds::mids2datlist(impdat_mids)
  dat_mice <- mice::complete(impdat_mids, "all")
  res_est <- evaluate_coxest_PH_multiple(impdat_mids, truelogHR)
  res_pred <- purrr::map(impdat_list, evaluate_coxperf_PH) %>% dplyr::bind_rows() %>% colMeans()
  cindex <-purrr::map(impdat_list, extract_cindex_PH) %>% dplyr::bind_rows() %>% colMeans()
  res_PH <- cbind(N, res_est, data.frame(t(res_pred)), cindex,  method, miss, mech, cens, setting = "PH", run_id = run_id, row.names = NULL)
  return(res_PH)
}

simcox_PH_multiple2 <- purrr::possibly(.f = simcox_PH_multiple, otherwise = NULL)
out_PH_multiple <- pbapply::pblapply(1:R, function(run_id) {
  purrr::map_df(grid_list, ~ simcox_PH_multiple2(
    N = .x$N, 
    mech = .x$mech, 
    miss = .x$miss, 
    method = .x$method, 
    cens = .x$cens, 
    run_id = run_id  # pass the loop variable as run_id
  ))
})
# Combining the results into a single data frame
res_PH_multiple <- do.call(rbind, out_PH_multiple)




## bind PH results -----
res_PH <- rbind(res_PH_single, res_PH_complete, res_PH_multiple)

write.csv(res_PH, "res_PH.csv")

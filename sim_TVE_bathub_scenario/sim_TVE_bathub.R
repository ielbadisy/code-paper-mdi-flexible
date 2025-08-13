#====================[ SETUP ]====================
rm(list = ls())
source("utils_TVE_bathub.R")

# Load required packages
pacman::p_load(
  dplyr, ggplot2, miceadds, rstpm2, stringr, biostat3, imputer, parallel, gower, 
  fastDummies, mice, survival, mixgb, survAUC, purrr, pbapply, tictoc, missCforest
)

#====================[ PARAMETERS ]====================
R <- 3
N <- 100
nCores <- 34
miss <- 0.3
cens <- 0.3
mech <- c("MNAR", "MAR", "MCAR")
shape <- "bathub"
set.seed(15102024)

#====================[ METHODS ]====================
valid_scenarios <- list(
  complete_data = "complete_data",
  complete = "complete",
  single = c("naive", "knn", "cart", "hotdeck", "famd", "missranger", "missforest", "misscforest"),
  multiple = c("mice", "micecart", "micerf")
)

grid_TVE <- do.call(rbind, lapply(names(valid_scenarios), function(scenario) {
  expand.grid(
    N = N,
    mech = mech,
    miss = miss,
    cens = cens,
    shape = shape,
    scenario = scenario,
    method = valid_scenarios[[scenario]],
    stringsAsFactors = FALSE
  )
}))

#====================[ MAIN SIMULATION FUNCTION ]====================
simcox_TVE <- function(N, mech, miss, method, cens, shape, run_id, scenario, ampdat) {
  dat <- genTVE(n = N, maxt = 5, cens = cens, shape = shape)
  res_obs <- obs_hr_x3(dat)
  res_est <- NULL
  
  if (scenario == "complete_data") {
    res_est <- evaluate_coxest_TVE(dat)
    
  } else if (scenario == "complete") {
    impdat <- suppressWarnings(imputer::imputer(ampdat, method = method))
    res_obs <- obs_hr_x3(impdat)
    res_est <- evaluate_coxest_TVE(impdat)
    
  } else if (scenario == "single") {
    impdat <- suppressWarnings(imputer::imputer(ampdat, method = method))
    res_est <- evaluate_coxest_TVE(impdat)
    
  } else if (scenario == "multiple") {
    impdat_mids <- suppressWarnings(imputer::imputer(ampdat, method = method))
    dat_mice <- complete(impdat_mids, "all")
    mod_imputation <- lapply(dat_mice, evaluate_coxest_TVE)
    
    mod_array <- simplify2array(mod_imputation)
    res_est <- as.data.frame(apply(mod_array, c(1, 2), mean))
    
    # Rubin's rules for CI pooling
    se_estimates <- (mod_array[, "hr_x3_upper", ] - mod_array[, "hr_x3_lower", ]) / (2 * 1.96)
    within_var <- apply(se_estimates^2, 1, mean)
    between_var <- apply(mod_array[, "hr_x3", ], 1, var)
    total_var <- within_var + (1 + 1/length(mod_imputation)) * between_var
    
    res_est$obs_hr_x3_lower <- res_est$hr_x3 - 1.96 * sqrt(total_var)
    res_est$obs_hr_x3_upper <- res_est$hr_x3 + 1.96 * sqrt(total_var)
  }
  
  # Final result assembly
  res_TVE <- data.frame(
    N, res_obs, res_est, method, miss, mech, cens, shape, setting = "TVE", run_id = run_id, row.names = NULL
  )
  
  # Empirical coverage for all cases
  # Empirical coverage
  res_TVE$coverage_hr_x3 <- if_else(
    res_TVE$method == "complete",
    data.table::between(res_TVE$obs_hr_x3, res_TVE$hr_x3_lower, res_TVE$hr_x3_upper),
    data.table::between(res_TVE$hr_x3, res_TVE$obs_hr_x3_lower, res_TVE$obs_hr_x3_upper)
  )
  
  return(res_TVE)
}

#====================[ RUN SIMULATIONS ]====================
grid_list <- pbapply::pblapply(seq_len(nrow(grid_TVE)), function(i) as.list(grid_TVE[i, ]))
simcox_TVE_safe <- purrr::safely(.f = simcox_TVE)

out_TVE <- pbapply::pblapply(1:R, function(run_id) {
  dat <- genTVE(n = N, maxt = 5, cens = cens, shape = shape)
  
  purrr::map_df(grid_list, function(param) {
    ampdat <- generate_NA_TVE(dat, mech = param$mech, pmiss = miss)
    
    result <- simcox_TVE_safe(
      N = param$N,
      mech = param$mech,
      miss = param$miss,
      method = param$method,
      shape = param$shape,
      cens = param$cens,
      run_id = run_id,
      scenario = param$scenario,
      ampdat = ampdat
    )
    
    if (!is.null(result$result)) result$result else NULL
  })
})

#====================[ EXPORT RESULTS ]====================
res_TVE <- do.call(rbind, out_TVE)
write.csv(res_TVE, "res_TVE_bathub.csv", row.names = FALSE)

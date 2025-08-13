#====================[ SETUP ]====================
rm(list = ls())
source("utils_NLE.R")

# Load packages
pacman::p_load(
  dplyr, ggplot2, miceadds, rstpm2, stringr, biostat3, imputer, parallel, gower, 
  fastDummies, mice, survival, mixgb, survAUC, purrr, pbapply, tictoc, missCforest
)

#====================[ PARAMETERS ]====================
R <- 500
N <- 1000
nCores <- 32
miss <- 0.3
cens <- 0.3
mech <- c("MNAR", "MAR", "MCAR")
set.seed(15102024)

#====================[ METHODS ]====================
valid_scenarios <- list(
  complete_data = "complete_data",
  complete = "complete",
  single = c("naive", "knn", "cart", "hotdeck", "famd", "missranger", "missforest", "misscforest"),
  multiple = c("mice", "micecart", "micerf")
)

grid_NLE <- do.call(rbind, lapply(names(valid_scenarios), function(scenario) {
  expand.grid(
    N = N,
    mech = mech,
    miss = miss,
    cens = cens,
    scenario = scenario,
    method = valid_scenarios[[scenario]],
    stringsAsFactors = FALSE
  )
}))

#====================[ CORE FUNCTION ]====================
simcox_NLE <- function(N, mech, miss, method, cens, run_id, scenario, ampdat) {
  dat <- genNLE(n = N, maxt = 5, cens = cens)
  res_obs <- obs_hr_x2(dat)
  res_est <- NULL
  
  if (scenario == "complete_data") {
    res_est <- evaluate_coxest_NLE(dat)
    
  } else if (scenario == "complete") {
    impdat <- suppressWarnings(imputer::imputer(ampdat, method = method))
    res_obs <- obs_hr_x2(impdat)
    res_est <- evaluate_coxest_NLE(impdat)
    
  } else if (scenario == "single") {
    impdat <- suppressWarnings(imputer::imputer(ampdat, method = method))
    res_est <- evaluate_coxest_NLE(impdat)
    
  } else if (scenario == "multiple") {
    impdat_mids <- suppressWarnings(imputer::imputer(ampdat, method = method))
    dat_mice <- complete(impdat_mids, "all")
    mod_imputation <- lapply(dat_mice, evaluate_coxest_NLE)
    
    mod_array <- simplify2array(mod_imputation)
    res_est <- as.data.frame(apply(mod_array, c(1, 2), mean))
    
    # Rubin's Rules pooling of CI
    se_estimates <- (mod_array[, "hr_x2_upper", ] - mod_array[, "hr_x2_lower", ]) / (2 * 1.96)
    within_var <- apply(se_estimates^2, 1, mean)
    between_var <- apply(mod_array[, "hr_x2", ], 1, var)
    total_var <- within_var + (1 + 1/length(mod_imputation)) * between_var
    
    res_est$obs_hr_x2_lower <- res_est$hr_x2 - 1.96 * sqrt(total_var)
    res_est$obs_hr_x2_upper <- res_est$hr_x2 + 1.96 * sqrt(total_var)
  }
  
  # Build results
  res_NLE <- data.frame(
    N, res_obs, res_est, method, miss, mech, cens,
    setting = "NLE", run_id = run_id, row.names = NULL
  )
  
  #====================[ EMPIRICAL COVERAGE ]====================
  res_NLE$coverage_hr_x2 <- if_else(
    res_NLE$method == "complete",
    data.table::between(res_NLE$obs_hr_x2, res_NLE$hr_x2_lower, res_NLE$hr_x2_upper),
    data.table::between(res_NLE$hr_x2, res_NLE$obs_hr_x2_lower, res_NLE$obs_hr_x2_upper)
  )
  
  return(res_NLE)
}

#====================[ PARALLEL EXECUTION ]====================
grid_list <- pbapply::pblapply(seq_len(nrow(grid_NLE)), function(i) as.list(grid_NLE[i, ]))
simcox_NLE_safe <- purrr::safely(.f = simcox_NLE)

out_NLE <- pbapply::pblapply(1:R, function(run_id) {
  dat <- genNLE(n = N, maxt = 5, cens = cens)
  dat <- imputer::imputer(dat, "complete")  # Remove missing values before generating new ones
  
  purrr::map_df(grid_list, function(param) {
    ampdat <- generate_NA_NLE(dat, mech = param$mech, pmiss = miss)
    result <- simcox_NLE_safe(
      N = param$N,
      mech = param$mech,
      miss = param$miss,
      method = param$method,
      cens = param$cens,
      run_id = run_id,
      scenario = param$scenario,
      ampdat = ampdat
    )
    if (!is.null(result$result)) result$result else NULL
  })
})

#====================[ OUTPUT HANDLING ]====================


res_NLE <- data.table::rbindlist(out_NLE[!sapply(out_NLE, is.null)], fill = TRUE)

# Save results
write.csv(res_NLE, "res_NLE.csv", row.names = FALSE)

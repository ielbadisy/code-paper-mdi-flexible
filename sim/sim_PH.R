# Scenario 1 (proportional hazards): 500 replicates x 3 mechanisms, N = 1000.
# In each replicate, one complete dataset is generated and amputed once per
# mechanism; all methods are applied to the same amputed dataset.
# Output: results/sim/res_PH.rds (one row per replicate x method x covariate).
# Used by analysis/res_PH.R (Tables 2-3, Appendix A).
source("R/labels.R")
source("R/utils_PH.R")
library(mcstatsim)

R_REPS <- 500
N_OBS <- 1000
MISS <- 0.3
N_CORES <- max(1, parallel::detectCores() - 1)
METHODS <- c("complete", "naive", "hotdeck", "knn", "cart", "famd", "missforest",
             "missranger", "misscforest", "mice", "micecart", "micerf")

#' One replicate of Scenario 1 for one missingness mechanism
#'
#' @param mech "MCAR", "MAR" or "MNAR"
#' @return data.frame of performance rows for all methods
simcox_PH <- function(mech) {
  dat <- genPH(n = N_OBS)
  ampdat <- generate_NA_PH(dat, mech = mech, pmiss = MISS)
  out <- lapply(METHODS, function(method) {
    res <- tryCatch({
      imp <- suppressWarnings(imputer::imputer(ampdat, method = method))
      if (method %in% MULTI_METHODS) evaluate_coxest_PH_multiple(imp)
      else evaluate_coxest_PH(imp)
    }, error = function(e) NULL)
    if (!is.null(res)) res$method <- method
    res
  })
  do.call(rbind, out)
}

grid_PH <- data.frame(mech = c("MNAR", "MAR", "MCAR"))
res_PH <- runsim(R_REPS, grid_PH, simcox_PH, num_cores = N_CORES, seed = 15102024,
                 show_progress = FALSE, on_error = "warn")
print(failure_summary(res_PH))
saveRDS(res_PH, results_path("sim", "res_PH.rds"))

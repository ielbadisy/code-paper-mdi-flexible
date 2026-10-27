#====================[ SETUP ]====================
# Scenario 5 (non-linear effect of x2). Runs all methods on the same amputed
# dataset within each replicate and saves per-replicate curves to res_NLE.rds.
# Requires mcstatsim >= 0.9 (remotes::install_github("ielbadisy/mcstatsim")).
rm(list = ls())
source("utils_NLE.R")
library(mcstatsim)

#====================[ PARAMETERS ]====================
R      <- 500
N      <- 1000
miss   <- 0.3
nCores <- 12

METHODS <- c("complete_data", "complete", "naive", "hotdeck", "knn", "cart", "famd",
             "missforest", "missranger", "misscforest", "mice", "micecart", "micerf")
MULTI   <- c("mice", "micecart", "micerf")

#====================[ ONE REPLICATE ]====================
simcox_NLE <- function(mech) {
  dat    <- genNLE(n = N)
  ampdat <- generate_NA_NLE(dat, mech = mech, pmiss = miss)
  out <- lapply(METHODS, function(method) {
    res <- tryCatch({
      if (method == "complete_data") {
        ci_single(fit_NLE(dat))
      } else if (method %in% MULTI) {
        imp <- suppressWarnings(imputer::imputer(ampdat, method = method))
        ci_rubin(lapply(mice::complete(imp, "all"), fit_NLE))
      } else {
        ci_single(fit_NLE(suppressWarnings(imputer::imputer(ampdat, method = method))))
      }
    }, error = function(e) data.frame(x2 = NA, est = NA, lower = NA, upper = NA))
    res$method <- method
    res
  })
  res <- do.call(rbind, out)
  res$truth <- g_x2(res$x2) - g_x2(X_REF)
  res$prop_missing <- mean(is.na(ampdat$x2))
  res
}

#====================[ RUN ]====================
grid_NLE <- data.frame(mech = c("MNAR", "MAR", "MCAR"))
res_NLE  <- runsim(R, grid_NLE, simcox_NLE, num_cores = nCores, seed = 15102024,
                   show_progress = FALSE, on_error = "warn",
                   save_path = "res_NLE_checkpoint.rds", checkpoint_every = 50L)
print(failure_summary(res_NLE))
saveRDS(res_NLE, "res_NLE.rds")

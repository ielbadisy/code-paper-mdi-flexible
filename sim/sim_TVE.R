# Scenarios 2-4: time-varying effect of x3 (bathtub, decreasing, increasing).
# Usage: Rscript sim/sim_TVE.R <shape> [n_replicates]
# Writes results/sim/res_TVE_<shape>.rds (one row per replicate x method x mechanism x time point).
#
# Replicates run serially after a single set.seed(15102024), in the order of the original
# run (about 30 seconds per replicate, about 4 hours per shape).

source("R/labels.R")
source("R/utils_TVE.R")
library(imputer)

args <- commandArgs(trailingOnly = TRUE)
SHAPE <- if (length(args) >= 1) args[1] else "increasing"
R_REPS <- if (length(args) >= 2) as.integer(args[2]) else 500
stopifnot(SHAPE %in% names(TVE_BETA))
N_OBS <- 1000
MISS <- 0.3
DF <- TVE_DF[[SHAPE]]

# Order of the conditions within a replicate (mechanism varies fastest)
METHODS <- list(
  complete_data = "complete_data",
  complete = "complete",
  single = c("naive", "knn", "cart", "hotdeck", "famd", "missranger", "missforest",
             "misscforest"),
  multiple = c("mice", "micecart", "micerf")
)
grid <- do.call(rbind, lapply(names(METHODS), function(sc) {
  expand.grid(mech = c("MNAR", "MAR", "MCAR"), scenario = sc, method = METHODS[[sc]],
              stringsAsFactors = FALSE)
}))

#' Estimate the hazard ratio curve of x3 for one condition
#'
#' @param ampdat data of the replicate with missing values in x3
#' @param scenario "complete_data", "complete", "single" or "multiple"
#' @param method imputation method passed to imputer::imputer()
simcox_TVE <- function(ampdat, scenario, method) {
  # The full-data reference (FullData) is fitted on a fresh complete dataset from the same
  # model. It is drawn for every condition, as in the original run, so that the
  # random-number stream is unchanged.
  fulldat <- genTVE(n = N_OBS, maxt = 5, shape = SHAPE)
  switch(scenario,
    complete_data = fit_TVE(fulldat, DF),
    complete = ,
    single = fit_TVE(suppressWarnings(imputer(ampdat, method = method)), DF),
    multiple = {
      mids <- suppressWarnings(imputer(ampdat, method = method))
      pool_TVE(lapply(mice::complete(mids, "all"), fit_TVE, df = DF))
    }
  )
}

set.seed(15102024)
out <- vector("list", R_REPS)
for (r in seq_len(R_REPS)) {
  dat <- genTVE(n = N_OBS, maxt = 5, shape = SHAPE)
  out[[r]] <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    ampdat <- generate_NA_TVE(dat, g$mech, MISS)
    est <- tryCatch(simcox_TVE(ampdat, g$scenario, g$method), error = function(e) NULL)
    if (is.null(est)) return(NULL)
    data.frame(rep = r, shape = SHAPE, mech = g$mech, method = g$method, est)
  }))
  message(sprintf("%s: replicate %d/%d done", SHAPE, r, R_REPS))
}
res <- do.call(rbind, out)
saveRDS(res, results_path("sim", sprintf("res_TVE_%s.rds", SHAPE)))

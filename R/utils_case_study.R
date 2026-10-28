# Helpers shared by the two colorectal cancer (CRC) case-study scripts.

#' Imputation methods compared in the case study (FullData is the reference)
CASE_METHODS <- c("complete", "naive", "hotdeck", "mice", "cart", "famd", "knn",
                  "misscforest", "missforest", "missranger", "micecart", "micerf")

#' Read and recode the CRC cohort (complete records of the analysis variables)
read_CRC <- function(path = file.path("data", "CRC.csv")) {
  dat <- imputer::imputer(read.csv(path), "complete")
  dat <- dat[, c("time", "event", "Age", "Surgery", "MSI.MSS.status", "Sex", "Stage")]
  dat$MSI.MSS.status <- as.numeric(dat$MSI.MSS.status == "Yes")
  dat$Stage <- as.numeric(as.factor(dat$Stage))
  dat$Surgery <- as.numeric(dat$Surgery == "Yes")
  dat$Sex <- as.numeric(dat$Sex == "Male")
  dat
}

#' Introduce missing values in one variable with mice::ampute
#'
#' @param dat CRC data
#' @param var name of the variable made incomplete
#' @param mech missingness mechanism
#' @param pmiss proportion of incomplete rows
ampute_CRC <- function(dat, var, mech = "MNAR", pmiss = 0.3) {
  patterns <- as.integer(names(dat) != var)
  mice::ampute(dat, mech = mech, prop = pmiss, patterns = patterns)$amp
}

#' Apply one imputation method to the amputed data and evaluate the model
#'
#' Single imputation gives one curve. Multiple imputation gives the average
#' of the curves over the imputed datasets.
#'
#' @param ampdat amputed data (the same for all methods)
#' @param method imputer method code
#' @param evaluate function(data) returning a data.frame with a numeric curve
impute_and_evaluate <- function(ampdat, method, evaluate) {
  if (method %in% MULTI_METHODS) {
    curves <- lapply(mice::complete(imputer::imputer(ampdat, method), "all"), evaluate)
    res <- Reduce(`+`, curves) / length(curves)
  } else {
    res <- evaluate(imputer::imputer(ampdat, method))
  }
  data.frame(res, method = method)
}

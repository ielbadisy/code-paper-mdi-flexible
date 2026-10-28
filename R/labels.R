# Method names used in the code and the labels used in the paper figures and tables.

#' Paper labels of the imputation methods, in the order of the figure panels
METHOD_LABELS <- c(
  complete_data = "FullData", complete = "CompleteCases", naive = "Median/mode",
  hotdeck = "Hotdeck", mice = "MICE", cart = "CART", famd = "FAMD", knn = "KNN",
  misscforest = "missCforest", missforest = "missForest", missranger = "missRanger",
  micecart = "miceCART", micerf = "miceRF"
)

#' Multiple imputation methods (pooled with Rubin's rules)
MULTI_METHODS <- c("mice", "micecart", "micerf")

#' Convert method codes to an ordered factor of paper labels
#'
#' @param x character vector of method codes (names of METHOD_LABELS)
#' @param full_data keep the full-data reference level ("FullData")
#' @return factor with levels in the order of METHOD_LABELS
method_factor <- function(x, full_data = FALSE) {
  lev <- if (full_data) METHOD_LABELS else METHOD_LABELS[-1]
  factor(METHOD_LABELS[as.character(x)], levels = lev)
}

#' Save a figure into results/figures
#'
#' @param plot ggplot object
#' @param name file name (without folder), e.g. "Figure1.png"
#' @param width,height size in inches
save_figure <- function(plot, name, width = 8, height = 6) {
  dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path("results", "figures", name), plot,
                  width = width, height = height, dpi = 300)
}

#' Path of a file in results/<sub>, creating the folder if needed
#'
#' @param sub sub-folder of results ("sim", "tables", "figures")
#' @param name file name
results_path <- function(sub, name) {
  dir.create(file.path("results", sub), recursive = TRUE, showWarnings = FALSE)
  file.path("results", sub, name)
}

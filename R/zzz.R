#' @import data.table
#' @importFrom magrittr "%>%"
NULL

rBiasCorrection_default_options <- list( # nolint
  rBiasCorrection.nls_implementation = "GN.paper" # nolint
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(rBiasCorrection_default_options) %in% names(op))
  if (any(toset)) options(rBiasCorrection_default_options[toset])
  invisible()
}

# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2020 Lorenz Kapsner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title better_model helper function
#'
#' @description Internal function to select the better model between
#'   hyperbolic regression and cubic regression.
#'
#' @inheritParams createbarerrorplots
#' @param statstable_post_hyperbolic A data.table object, containing
#'   the output of \code{statisticsList_()} of the calculated regression
#'   parameters form the calibration data corrected with hyperbolic regression.
#' @param statstable_post_cubic A data.table object, containing the output
#'   of \code{statisticsList_()} of the calculated regression parameters form
#'   the calibration data corrected with cubic regression.
#' @inheritParams biascorrection
#'
#' @return The function returns a data.table with 4 columns, the last column
#'   being named 'better_model', which indicates in a binary manner,
#'   if the hyperbolic model (better_model = 0) or the cubic model
#'   (better_model = 1) result in a 'better' `SSE` or `RelError` respectively.
#'
#' @examples
#' \donttest{
#' # define list object to save all data
#' rv <- list()
#' rv$minmax <- TRUE
#' rv$selection_method <- "RelError"
#' rv$sample_locus_name <- "Test"
#' rv$seed <- 1234
#'
#' # define logfilename
#' logfilename <- paste0(tempdir(), "/log.txt")
#'
#' # import experimental file
#' exp_type_1 <- rBiasCorrection::example.data_experimental
#' rv$fileimport_experimental <- exp_type_1$dat
#'
#' # import calibration file
#' cal_type_1 <- rBiasCorrection::example.data_calibration
#' rv$fileimport_calibration <- cal_type_1$dat
#' rv$vec_cal <- cal_type_1$vec_cal
#'
#'
#' # perform regression
#' regression_results <- regression_utility(
#'   rv$fileimport_calibration,
#'   "Testlocus",
#'   locus_id = NULL,
#'   rv = rv,
#'   mode = NULL,
#'   logfilename,
#'   minmax = rv$minmax,
#'   seed = rv$seed
#' )
#'
#' # extract regression results
#' rv$result_list <- regression_results$result_list
#'
#' # get regression statistics
#' rv$reg_stats <- statistics_list(
#'   rv$result_list,
#'   minmax = rv$minmax
#' )
#'
#' # select the better model based on the sum of squared errrors ("SSE")
#' rv$choices_list <- better_model(
#'   statstable_pre = rv$reg_stats,
#'   selection_method = "SSE"
#' )
#' }
#'
#' @export
#'
better_model <- function(statstable_pre,
                         statstable_post_hyperbolic = NULL,
                         statstable_post_cubic = NULL,
                         selection_method = "SSE") {

  stopifnot(
    is.character(selection_method),
    selection_method %in% c("SSE", "RelError"),
    data.table::is.data.table(statstable_pre),
    data.table::is.data.table(statstable_post_hyperbolic) ||
      is.null(statstable_post_hyperbolic),
    data.table::is.data.table(statstable_post_cubic) ||
      is.null(statstable_post_cubic)
  )

  if (selection_method == "SSE") {
    outdat <- statstable_pre[, c("Name",
                                "SSE_hyperbolic",
                                "SSE_cubic"), with = F]
    # mark the better model: 1 = cubic, 0 = hyperbolic
    outdat[, ("better_model") := ifelse(
      get("SSE_cubic") <= get("SSE_hyperbolic"),
      1,
      0)]
  } else if (selection_method == "RelError") {
    x_dat <- statstable_post_hyperbolic[, c("Name", "relative_error"), with = F]
    colnames(x_dat) <- c("Name", "relative_error_h")
    y_dat <- statstable_post_cubic[, c("Name", "relative_error"), with = F]
    outdat <- merge(x = x_dat,
                    y = y_dat,
                    by = "Name",
                    all = T,
                    suffixes = c("", "_c"))
    colnames(outdat) <- c("Name", "relative_error_h", "relative_error_c")
    outdat[, ("better_model") := ifelse(
      get("relative_error_c") <= get("relative_error_h"),
      1,
      0)]
  }
  return(outdat)
}

# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2025 Lorenz Kapsner
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

#' @title statistics_list helper function
#'
#' @description Internal function that converts the results_list (output of
#'   \code{regression_utility()}) into a data.table object.
#'
#' @param resultlist A list object. The results_list output of
#'   \code{regression_utility()}.
#' @inheritParams biascorrection
#'
#' @return The function takes the `resultslist` and converts it to a
#'   statistics list, which is basically a `data.table` with the results
#'   of the hyperbolic and the cubic regression.
#'
#' @examples
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
#' @export
#'
statistics_list <- function(resultlist,
                            minmax = FALSE) {
  if (isFALSE(minmax)) {
    outdat <- data.table::data.table(
      "Name" = character(),
      "relative_error" = numeric(),
      "SSE_hyperbolic" = numeric(),
      "R2_hyperbolic" = numeric(),
      "a_hyperbolic" = numeric(),
      "b_hyperbolic" = numeric(),
      "d_hyperbolic" = numeric(),
      "b1_hyperbolic" = numeric(),
      "s_hyperbolic" = numeric(),
      "###" = character(),
      "SSE_cubic" = numeric(),
      "R2_cubic" = numeric(),
      "a_cubic" = numeric(),
      "b_cubic" = numeric(),
      "c_cubic" = numeric(),
      "d_cubic" = numeric()
    )

    out_names <- c(
      "Name",
      "relative_error",
      "SSE_hyperbolic",
      "R2_hyperbolic",
      "a_hyperbolic",
      "b_hyperbolic",
      "d_hyperbolic",
      "b1_hyperbolic",
      "s_hyperbolic",
      "SSE_cubic",
      "R2_cubic",
      "a_cubic",
      "b_cubic",
      "c_cubic",
      "d_cubic"
    )

  } else if (isTRUE(minmax)) {
    outdat <- data.table::data.table(
      "Name" = character(),
      "relative_error" = numeric(),
      "SSE_hyperbolic" = numeric(),
      "R2_hyperbolic" = numeric(),
      "b_hyperbolic" = numeric(),
      "###" = character(),
      "SSE_cubic" = numeric(),
      "R2_cubic" = numeric(),
      "a_cubic" = numeric(),
      "b_cubic" = numeric(),
      "####" = character(),
      "y0" = numeric(),
      "y1" = numeric(),
      "m0" = numeric(),
      "m1" = numeric()
    )

    out_names <- c(
      "Name",
      "relative_error",
      "SSE_hyperbolic",
      "R2_hyperbolic",
      "b_hyperbolic",
      "SSE_cubic",
      "R2_cubic",
      "a_cubic",
      "b_cubic",
      "y0",
      "y1",
      "m0",
      "m1"
    )
  }

  for (i in names(resultlist)) {

    if (isFALSE(minmax)) {
      out_list <- list(
        i,
        resultlist[[i]][["relative_error"]],
        resultlist[[i]][["SSE_hyper"]],
        resultlist[[i]][["Coef_hyper"]][["R2"]],
        resultlist[[i]][["Coef_hyper"]][["a"]],
        resultlist[[i]][["Coef_hyper"]][["b"]],
        resultlist[[i]][["Coef_hyper"]][["d"]],
        resultlist[[i]][["Coef_hyper"]][["b1"]],
        resultlist[[i]][["Coef_hyper"]][["s"]],
        resultlist[[i]][["SSE_cubic"]],
        resultlist[[i]][["Coef_cubic"]][["R2"]],
        resultlist[[i]][["Coef_cubic"]][["a"]],
        resultlist[[i]][["Coef_cubic"]][["b"]],
        resultlist[[i]][["Coef_cubic"]][["c"]],
        resultlist[[i]][["Coef_cubic"]][["d"]]
      )
    } else if (isTRUE(minmax)) {
      out_list <- list(
        i,
        resultlist[[i]][["relative_error"]],
        resultlist[[i]][["SSE_hyper"]],
        resultlist[[i]][["Coef_hyper"]][["R2"]],
        resultlist[[i]][["Coef_hyper"]][["b"]],
        resultlist[[i]][["SSE_cubic"]],
        resultlist[[i]][["Coef_cubic"]][["R2"]],
        resultlist[[i]][["Coef_cubic"]][["a"]],
        resultlist[[i]][["Coef_cubic"]][["b"]],
        resultlist[[i]][["Coef_hyper"]][["y0"]],
        resultlist[[i]][["Coef_hyper"]][["y1"]],
        resultlist[[i]][["Coef_hyper"]][["m0"]],
        resultlist[[i]][["Coef_hyper"]][["m1"]]
      )
    }

    names(out_list) <- out_names

    out_list <- sapply(
      X = out_list,
      FUN = function(x) {
        if (!is.character(x)) {
          x <- as.numeric(as.character(x))
        }
        return(x)
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    outdat <- data.table::rbindlist(
      l = list(
        outdat,
        out_list
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }
  #% mark the better model: 1 = cubic, 0 = hyperbolic
  #% outdat[, ("better_model") := ifelse(
  #%   get("SSE_cubic") <= get("SSE_hyperbolic"),
  #%   1,
  #%   0
  #% )]
  return(outdat)
}

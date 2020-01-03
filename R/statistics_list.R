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

#' @title statistics_list helper function
#'
#' @description Internal function that converts the results_list (output of
#'   \code{regressionUtility_()}) into a data.table object.
#'
#' @param resultlist A list object. The results_list output of
#'   \code{regressionUtility_()}.
#' @inheritParams biascorrection
#'
#' @export
#'
statistics_list <- function(resultlist,
                            minmax = FALSE) {
  if (isFALSE(minmax)) {
    outdat <- data.table::data.table("Name" = names(resultlist),
                                     "relative_error" = NA,
                                     "SSE_hyperbolic" = NA,
                                     "R2_hyperbolic" = NA,
                                     "a_hyperbolic" = NA,
                                     "b_hyperbolic" = NA,
                                     "d_hyperbolic" = NA,
                                     "b1_hyperbolic" = NA,
                                     "s_hyperbolic" = NA,
                                     "###" = NA,
                                     "SSE_cubic" = NA,
                                     "R2_cubic" = NA,
                                     "a_cubic" = NA,
                                     "b_cubic" = NA,
                                     "c_cubic" = NA,
                                     "d_cubic" = NA)

    outdat[, ("Name") := names(resultlist)]

    vec <- names(outdat)[-1]
    outdat[, (vec) := lapply(.SD, function(x) {
      as.numeric(as.character(x))
      }), .SDcols = vec]

    for (i in names(resultlist)) {
      out_names <- c("relative_error",
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
                     "d_cubic")
      out_list <- list(resultlist[[i]][["relative_error"]],
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
                       resultlist[[i]][["Coef_cubic"]][["d"]])
      outdat[get("Name") == i, (out_names) := out_list]
    }
  } else if (isTRUE(minmax)) {
    outdat <- data.table::data.table("Name" = names(resultlist),
                                     "relative_error" = NA,
                                     "SSE_hyperbolic" = NA,
                                     "R2_hyperbolic" = NA,
                                     "b_hyperbolic" = NA,
                                     "###" = NA,
                                     "SSE_cubic" = NA,
                                     "R2_cubic" = NA,
                                     "a_cubic" = NA,
                                     "b_cubic" = NA,
                                     "####" = NA,
                                     "y0" = NA,
                                     "y1" = NA,
                                     "m0" = NA,
                                     "m1" = NA)

    outdat[, ("Name") := names(resultlist)]

    vec <- names(outdat)[-1]
    outdat[, (vec) := lapply(.SD, function(x) {
      as.numeric(as.character(x))
      }), .SDcols = vec]

    for (i in names(resultlist)) {
      out_names <- c("relative_error",
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
                     "m1")
      out_list <- list(resultlist[[i]][["relative_error"]],
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
                       resultlist[[i]][["Coef_hyper"]][["m1"]])
      outdat[get("Name") == i, (out_names) := out_list]
    }
  }
  #% mark the better model: 1 = cubic, 0 = hyperbolic
  #% outdat[, ("better_model") := ifelse(
  #%   get("SSE_cubic") <= get("SSE_hyperbolic"),
  #%   1,
  #%   0
  #% )]
  return(outdat)
}

# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2022 Lorenz Kapsner
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

#' @title aggregated_input helper function
#'
#' @description Internal function to present aggregated input data on which
#'   calculations are performed. This function does only have an effect, if
#'   repeated measurements are used for calibration a/o experimental data.
#'
#' @param vec_cal The vector containing the CpG columns (output of `clean_dt`
#'   with description = "calibration").
#'
#' @return A data.table in the long format with aggregated means for each CpG
#'   site of each sample and the corresponding standard deviation.
#'
#' @inheritParams clean_dt
#'
#' @examples
#' experimental <- rBiasCorrection::example.data_experimental
#' calibration <- rBiasCorrection::example.data_calibration
#'
#' vec_cal <- calibration$vec_cal
#'
#' experimental_aggregated <- aggregated_input(
#'   datatable = experimental$dat,
#'   description = "experimental",
#'   vec_cal = vec_cal,
#'   type = 1
#' )
#' dim(experimental_aggregated)
#' class(experimental_aggregated)
#'
#' calibration_aggregated <- aggregated_input(
#'   datatable = calibration$dat,
#'   description = "calibration",
#'   vec_cal = vec_cal
#' )
#' dim(calibration_aggregated)
#' class(calibration_aggregated)
#'
#' @export
#'
aggregated_input <- function(
    datatable,
    description = c("experimental", "calibration"),
    vec_cal,
    type = NULL
  ) {

  description <- match.arg(description)

  if (!is.null(type)) {
    type <- as.integer(type)
  }

  stopifnot(
    "`datatable` must be of class data.table" =
      data.table::is.data.table(datatable),
    is.vector(vec_cal),
    ifelse(
      test = !is.null(type),
      yes = description == "experimental" && type == 1L || type == 2L,
      no = TRUE
    )
  )

  outdat <- data.table::data.table()

  if (description == "experimental") {
    for (i in vec_cal) {
      append_df <- create_agg_df_exp(
        datatable = datatable,
        index = i,
        type = type
      )
      colnames(append_df)[2] <- "mean"
      # add CpG name
      append_df <- cbind(
        "CpG" = i,
        append_df
      )

      # concat to outdat
      outdat <- data.table::rbindlist(
        l = list(
          outdat,
          append_df
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
  } else if (description == "calibration") {
    for (i in vec_cal) {
      append_df <- create_agg_df(
        datatable = datatable,
        index = i
      )
      colnames(append_df)[2] <- "mean"
      # add CpG name
      append_df <- cbind(
        "CpG" = i,
        append_df
      )

      # concat to outdat
      outdat <- data.table::rbindlist(
        l = list(
          outdat,
          append_df
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
  }
  return(outdat)
}

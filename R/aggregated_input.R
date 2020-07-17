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
#'   site of each sample and the correponding standard deviation.
#'
#' @inheritParams clean_dt
#'
#' @export
#'
aggregated_input <- function(datatable, description, vec_cal, type = NULL) {

  stopifnot(
    data.table::is.data.table(datatable),
    is.character(description),
    description %in% c("experimental", "calibration"),
    is.vector(vec_cal),
    ifelse(
      description == "experimental",
      type == 1 || type == 2,
      is.null(type)
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

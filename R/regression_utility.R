# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2021 Lorenz Kapsner
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


#' @title regression_utility helper function
#'
#' @description Internal function to carry out the regression calculations.
#'
#' @param data A data.table object that contains the calibration data.
#' @inheritParams biascorrection
#' @param locus_id A character string. Default: NULL. ID of the respective
#'   locus (only used in type 2 correction).
#' @param rv A list object. A list that contains additional objects needed
#'   for the algorithms.
#' @param mode A character string. Default: NULL. Used to indicate "corrected"
#'   calibration data.
#'
#' @return The function performs the regression calculations and returns
#'   the results in a list.
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
#' length(regression_results)
#' class(regression_results)
#'
#' @export
#'
regression_utility <- function(data,
                               samplelocusname,
                               locus_id = NULL,
                               rv,
                               mode = NULL,
                               logfilename,
                               minmax,
                               seed = 1234) {

  if (!is.null(locus_id)) {
    write_log(
      message = paste0("### Starting with regression ",
                       "calculations ###\n\nLocus ID: ",
                       locus_id),
      logfilename = logfilename)
  } else {
    write_log(
      message = paste0("### Starting with regression calculations ###"),
      logfilename = logfilename)
  }

  regression_results <- regression_type1(
    datatable = data,
    vec_cal = rv$vec_cal,
    mode = mode,
    logfilename = logfilename,
    minmax = minmax,
    locus_id,
    locusname = rv$sample_locus_name,
    seed = seed
  )

  return(regression_results)
}

regression_type1 <- function(datatable,
                             vec_cal,
                             mode = NULL,
                             logfilename,
                             minmax,
                             locus_id = NULL,
                             locusname,
                             seed = 1234) {
  write_log(
    message = "Entered 'regression_type1'-Function",
    logfilename = logfilename
  )

  if (is.null(locus_id)) {
    plot_title <- "bquote(italic(.(locusname)))"
  } else {
    plot_title <- paste0(
      "bquote('Locus: ' ~ italic(.(locus_id)) ~",
      "' - Sample: ' ~ .(locusname))"
    )
  }

  # result_list
  result_list <- future.apply::future_sapply(
    X = vec_cal,
    FUN = function(i) {
      local({
        message <- paste0("# CpG-site: ", i)
        write_log(message = message,
                  logfilename = logfilename)
        df_agg <- create_agg_df(
          datatable = datatable,
          index = i
        )

        message(df_agg)
        write_log(
          message = paste("Logging df_agg:", i),
          logfilename = logfilename
        )
        write_log(message = df_agg, logfilename = logfilename)

        hr <- hyperbolic_regression(
          df_agg = df_agg,
          vec = i,
          logfilename = logfilename,
          minmax = minmax,
          seed = seed
        )

        cr <- cubic_regression(
          df_agg = df_agg,
          vec = i,
          logfilename = logfilename,
          minmax = minmax,
          seed = seed
        )
        # append result_list
        return(c(hr, cr))
      })
    },
    USE.NAMES = TRUE,
    simplify = FALSE,
    future.seed = TRUE
  )


  plot.listR <- future.apply::future_lapply(
    X = seq_len(length(vec_cal)),
    FUN = function(i) {
      local({

        df_agg <- create_agg_df(
          datatable = datatable,
          index = vec_cal[i]
        )
        if (is.null(mode)) {
          custom_ylab <- "methylation (%)\napparent after quantification"
        } else if (mode == "corrected") {
          custom_ylab <- "methylation (%)\nafter BiasCorrection"
        }

        lb1 <- c(paste0("  R\u00B2: \n  Hyperbolic = ",
                        round(result_list[[vec_cal[i]]]$Coef_hyper$R2, 2),
                        "\n  Cubic = ",
                        round(result_list[[vec_cal[i]]]$Coef_cubic$R2, 2)), "")


        gdat <- df_agg[
          , ("true_methylation") := as.numeric(
            as.character(
              get("true_methylation")
            )
          )
        ][
          , ("CpG") := as.numeric(
            as.character(
              get("CpG")
            )
          )
        ]

        if (is.null(mode)) {
          if ("sd" %in% colnames(gdat)) {
            gdat <- gdat[
              , ("sd") := as.numeric(
                as.character(
                  get("sd")
                )
              )
            ][
              , ("ymin") := get("CpG") - get("sd")
            ][
              , ("ymax") := get("CpG") + get("sd")
            ]
          }
        }

        p <- ggplot2::ggplot(data = gdat,
                             ggplot2::aes_string(
                               x = "true_methylation",
                               y = "CpG")
        ) +
          ggplot2::geom_point() +
          ggplot2::ylab(custom_ylab) +
          ggplot2::xlab("actual methylation (%)") +
          ggplot2::labs(
            title = eval(parse(text = plot_title)),
            subtitle = paste("CpG:", vec_cal[i])
          ) +
          ggplot2::geom_text(
            data = data.frame(),
            ggplot2::aes(x = -Inf,
                         y = c(max(gdat$true_methylation),
                               0.95 * max(gdat$true_methylation)),
                         hjust = 0, vjust = 1),
            label = lb1,
            parse = FALSE
          )

        if (is.null(mode)) {
          if ("ymin" %in% colnames(gdat) &&
              "ymax" %in% colnames(gdat)) {
            p <- p + ggplot2::geom_pointrange(
              ggplot2::aes_string(
                ymin = "ymin",
                ymax = "ymax"
              )
            )
          }
        }
        return(p)
      })
    },
    future.seed = TRUE
  )
  return(
    list("plot_list" = plot.listR,
         "result_list" = result_list)
  )
}

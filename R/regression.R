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
#' @param headless A logical (default: FALSE). Indicates, if the function is
#'   called from within a shiny app (default) or just from the console without
#'   a graphical user interface.
#' @param mode A character string. Default: NULL. Used to indicate "corrected"
#'   calibration data.
#'
#' @export
#'
regression_utility <- function(data,
                               samplelocusname,
                               locus_id = NULL,
                               rv,
                               mode = NULL,
                               headless = FALSE,
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


  # workaround to hide shiny-stuff, when going headless
  if (isFALSE(headless)) {
    # for plotting: basic idea and some code snippets from:
    # https://gist.github.com/wch/5436415/
    regression <- shiny::reactive({
      regression_type1(datatable = data,
                       vec_cal = rv$vec_cal,
                       mode = mode,
                       logfilename = logfilename,
                       minmax = minmax,
                       locus_id,
                       locusname = rv$sample_locus_name,
                       seed = seed)
    })

    shiny::withProgress(
      message = "Calculating calibration curves",
      value = 0,
      expr = {
        shiny::incProgress(1 / 1,
                           detail = "... working on calculations ...")
        # calculate results (if this is run here, j must be resetted)
        regression_results <- regression()
      }
    )
  } else {
    regression_results <- regression_type1(datatable = data,
                                           vec_cal = rv$vec_cal,
                                           mode = mode,
                                           logfilename = logfilename,
                                           minmax = minmax,
                                           locus_id,
                                           locusname = rv$sample_locus_name,
                                           seed = seed)
  }

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
    plot_title <- bquote(italic(.(locusname)))
  } else {
    plot_title <- bquote("Locus: " ~
                           italic(.(locus_id)) ~
                           " - Sample: " ~
                           .(locusname))
  }

  # result_list
  result_list <- list()

  plot.listR <- list()

  for (i in seq_len(length.out = length(vec_cal))) {
    message <- paste0("# CpG-site: ", vec_cal[i])
    write_log(message = message,
              logfilename = logfilename)
    df_agg <- stats::na.omit(
      create_agg_df(
        datatable = datatable,
        index = vec_cal[i]
      )
    )

    print(df_agg)
    write_log(
      message = paste("Logging df_agg:", vec_cal[i]),
      logfilename = logfilename
    )
    write_log(message = df_agg, logfilename = logfilename)

    hr <- hyperbolic_regression(
      df_agg = df_agg,
      vec = vec_cal[i],
      logfilename = logfilename,
      minmax = minmax,
      seed = seed
    )

    cr <- cubic_regression(
                                     df_agg = df_agg,
                                     vec = vec_cal[i],
                                     logfilename = logfilename,
                                     minmax = minmax,
                                     seed = seed
                                   )
    # append result_list
    result_list[[vec_cal[i]]] <- list(hr, cr)

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

    p <- ggplot2::ggplot(data = gdat,
                         ggplot2::aes_string(
                           x = "true_methylation",
                           y = "CpG")
    ) +
      ggplot2::geom_point() +
      ggplot2::ylab(custom_ylab) +
      ggplot2::xlab("actual methylation (%)") +
      ggplot2::labs(
        title = plot_title,
        subtitle = paste("CpG:", vec_cal[i])
      ) +
      ggplot2::geom_text(
        data = data.frame(),
        ggplot2::aes(x = -Inf,
                     y = c(max(gdat$true_methylation),
                           0.95 * max(gdat$true_methylation)),
                     hjust = 0, vjust = 1),
        label = lb1,
        parse = F
      )
    plot.listR[[i]] <- p
  }
  return(
    list("plot_list" = plot.listR,
         "result_list" = result_list)
  )
}

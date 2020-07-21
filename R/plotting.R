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

#' @title plotting_utility helper function
#'
#' @description Internal function to carry out the plotting of the
#'   calibrations curves.
#'
#' @param plotlist_reg A list object contating regression plots without
#'   regression curves (output of \code{regression_utility()}).
#' @inheritParams regression_utility
#' @inheritParams createbarerrorplots
#'
#' @return This function creates calibration plots and writes them to the
#'   local filesystem.
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
#' logfilename <- "log.txt"
#'
#' # define plotdir
#' rv$plotdir <- paste0(tempdir(), "/plots/")
#' dir.create(rv$plotdir)
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
#' # extract the plotlist
#' plotlist_reg <- regression_results$plot_list
#'
#' plotting_utility(
#'   data = rv$fileimport_calibration,
#'   plotlist_reg = plotlist_reg,
#'   type = 1,
#'   samplelocusname = rv$sample_locus_name,
#'   locus_id = NULL,
#'   rv = rv,
#'   mode = NULL,
#'   plotdir = rv$plotdir,
#'   logfilename = logfilename,
#'   minmax = rv$minmax,
#'   plot_height = 5,
#'   plot_width = 7.5,
#'   plot_textsize = 1
#' )
#'
#' @export
#'
# plotting utility
plotting_utility <- function(data,
                             plotlist_reg,
                             type,
                             samplelocusname,
                             locus_id = NULL,
                             rv,
                             mode = NULL,
                             plotdir,
                             logfilename,
                             minmax,
                             plot_height = 5,
                             plot_width = 7.5,
                             plot_textsize = 1) {

  if (!is.null(locus_id)) {
    write_log(
      message = paste0("### Starting with plotting ###\n\nLocus ID: ",
                       locus_id),
      logfilename = logfilename
    )
  } else {
    write_log(
      message = paste0("### Starting with plotting ###"),
      logfilename = logfilename
    )
  }

  # get number of CpG-sites
  length_vector <- length(rv$vec_cal)

  # get result_list
  if (is.null(mode)) {
    result_list <- rv$result_list
  } else if (mode == "corrected_h") {
    result_list <- rv$result_list_hyperbolic
  } else if (mode == "corrected_c") {
    result_list <- rv$result_list_cubic
  }

  vec_cal <- rv$vec_cal

  future.apply::future_Map(function(f) {
    local({
      plotname <- paste0(gsub("[[:punct:]]", "", vec_cal[f]))

      # filename-suffix
      fn_suffix <- ifelse(is.null(mode), "", paste0("_", mode))
      # message-suffix
      msg_suffix <- ifelse(is.null(mode), "", ifelse(
        mode == "corrected_h",
        "BiasCorrected (hyperbolic)",
        "BiasCorrected (cubic)")
      )

      # filname of temporary plot
      if (type == 1) {
        filename <- paste0(plotdir,
                           samplelocusname,
                           "_",
                           plotname,
                           fn_suffix,
                           ".png")
        plotmessage <- paste0("Creating ",
                              msg_suffix,
                              " plot No. ",
                              f)
      } else if (type == 2) {
        filename <- paste0(plotdir,
                           locus_id,
                           "-",
                           samplelocusname,
                           "_",
                           plotname,
                           fn_suffix,
                           ".png")
        plotmessage <- paste0("Locus ID: ",
                              locus_id,
                              " --> Creating ",
                              msg_suffix,
                              " plot No. ",
                              f)
      }

      write_log(
        message = paste(plotmessage, "- filename:", filename),
        logfilename = logfilename
      )

      # store plots to local temporary file
      create_plots(plotlist = plotlist_reg[[f]],
                   f = f,
                   vec_cal = vec_cal,
                   result_list = result_list,
                   filename = filename,
                   logfilename = logfilename,
                   mode = mode,
                   minmax = minmax,
                   plot_height = plot_height,
                   plot_width = plot_width,
                   plot_textsize = plot_textsize)
    })
  }, 1:length_vector)
}

create_plots <- function(plotlist,
                         f,
                         vec_cal,
                         result_list,
                         filename,
                         logfilename,
                         mode = NULL,
                         minmax,
                         plot_height = 5,
                         plot_width = 7.5,
                         plot_textsize = 16) {

  # hyperbolic parameters
  coef_h <- result_list[[vec_cal[f]]][["Coef_hyper"]]
  # cubic parameters
  coef_c <- result_list[[vec_cal[f]]][["Coef_cubic"]]

  if (isFALSE(minmax)) {

    # create messages
    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", a =", coef_h$a,
                  ", b =", coef_h$b,
                  ", d =", coef_h$d)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Cubic: Using a =", coef_c$a,
                  ", b =", coef_c$b,
                  ", c =", coef_c$c,
                  ", d =", coef_c$d)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    outplot <- calibration_plot(
      plotlist = plotlist,
      coef_hyper = coef_h,
      coef_cubic = coef_c,
      plot_textsize = plot_textsize,
      minmax = minmax
    )

  } else if (isTRUE(minmax)) {
    # create messages
    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", y0 =", coef_h$y0,
                  ", y1 =", coef_h$y1,
                  ", m0 =", coef_h$m0,
                  ", m1 =", coef_h$m1)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Cubic: Using a =", coef_c$a,
                  ", b =", coef_c$b,
                  ", y0 =", coef_c$y0,
                  ", y1 =", coef_c$y1,
                  ", m0 =", coef_c$m0,
                  ", m1 = ", coef_c$m1)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    outplot <- calibration_plot(
      plotlist = plotlist,
      coef_hyper = coef_h,
      coef_cubic = coef_c,
      plot_textsize = plot_textsize,
      minmax = minmax
    )
  }

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = plot_height,
    width = plot_width,
    dpi = 600
  )
}


#' @title createbarerrorplots helper function
#'
#' @description Internal function to create relative-error bar plots.
#'
#' @param statstable_pre A data.table object, containing the output of
#'   \code{statisticsList_()} of the calculated regression parameters (form
#'   the provided calibration data).
#' @param statstable_post A data.table object, containing the output of
#'   \code{statisticsList_()} of the calculated regression parameters form
#'   the corrected calibration data.
#' @param plot_height A integer value. The height (unit: inch) of the
#'   resulting plots (default: 5).
#' @param plot_width A integer value. The width (unit: inch) of the
#'   resulting plots (default: 7.5).
#' @param plot_textsize A integer value. The textsize of the
#'   resulting plots (default: 16).
#' @inheritParams regression_utility
#' @inheritParams clean_dt
#' @inheritParams on_start
#'
#' @return This function creates error bar-plots to visualize the realitve
#'   error before and after bias correction and writes these plots to the
#'   local filesystem.
#'
#' @examples
#' # define list object to save all data
#' rv <- list()
#' rv$minmax <- TRUE
#' rv$selection_method <- "RelError"
#' rv$sample_locus_name <- "Test"
#' rv$seed <- 1234
#'
#' # define plotdir
#' rv$plotdir <- paste0(tempdir(), "/plots/")
#' dir.create(rv$plotdir)
#'
#' # define logfilename
#' logfilename <- "log.txt"
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
#'   minmax = TRUE
#' )
#'
#' # select the better model based on the sum of squared errrors ("SSE")
#' rv$choices_list <- better_model(
#'   statstable_pre = rv$reg_stats,
#'   selection_method = "SSE"
#' )
#'
#' # correct calibration data (to show corrected calibration curves)
#' solved_eq_h <- solving_equations(datatable = rv$fileimport_calibration,
#'                                  regmethod = rv$choices_list,
#'                                  type = 1,
#'                                  rv = rv,
#'                                  mode = "corrected",
#'                                  logfilename = logfilename,
#'                                  minmax = rv$minmax)
#' rv$fileimport_cal_corrected_h <- solved_eq_h$results
#' colnames(rv$fileimport_cal_corrected_h) <- colnames(
#'   rv$fileimport_calibration
#' )
#'
#' # calculate new calibration curves from corrected calibration data
#' regression_results <- regression_utility(
#'   data = rv$fileimport_cal_corrected_h,
#'   samplelocusname = rv$sample_locus_name,
#'   rv = rv,
#'   mode = "corrected",
#'   logfilename = logfilename,
#'   minmax = rv$minmax,
#'   seed = rv$seed
#' )
#' rv$result_list_hyperbolic <- regression_results$result_list
#'
#'
#' # save regression statistics to reactive value
#' rv$reg_stats_corrected_h <- statistics_list(
#'   resultlist = rv$result_list_hyperbolic,
#'   minmax = rv$minmax
#' )
#'
#' createbarerrorplots(
#'   statstable_pre = rv$reg_stats,
#'   statstable_post = rv$reg_stats_corrected_h,
#'   rv = rv,
#'   type = 1,
#'   locus_id = NULL,
#'   plotdir = rv$plotdir,
#'   logfilename = logfilename,
#'   mode = "corrected_h",
#'   plot_height = 5,
#'   plot_width = 7.5,
#'   plot_textsize = 1
#' )
#'
#' @export
#'
createbarerrorplots <- function(statstable_pre,
                                statstable_post,
                                rv,
                                type,
                                locus_id = NULL,
                                plotdir,
                                logfilename,
                                mode = NULL,
                                plot_height = 5,
                                plot_width = 7.5,
                                plot_textsize = 16) {

  stats_pre <- statstable_pre[, c("Name", "relative_error"), with = F]
  stats_post <- statstable_post[, c("Name", "relative_error"), with = F]

  # calc ylim based on observed relative error
  max_err <- ifelse(max(stats_pre$relative_error) >
                      max(stats_post$relative_error),
                    max(stats_pre$relative_error),
                    max(stats_post$relative_error))

  ylim_max <- round_to_fifty(max_err)

  error_data <- merge(
    x = stats_post,
    y = stats_pre,
    by = "Name",
    sort = F,
    suffixes = c("", "_pre")
  )

  if (is.null(locus_id)) {
    locus <- rv$sample_locus_name
  } else {
    locus <- paste("Locus:",
                   locus_id,
                   "-",
                   "Sample:",
                   rv$sample_locus_name)
  }

  # Test if names are eqal
  if (identical(stats_pre[, get("Name")], stats_post[, get("Name")])) {

    vec_cal <- stats_pre[, get("Name")]
    length_vector <- length(vec_cal)

    sample_locus_name <- rv$sample_locus_name

    future.apply::future_Map(function(i) {
      local({
        plotname <- paste0(gsub("[[:punct:]]", "", vec_cal[i]))

        if (type == 1) {
          filename <- paste0(plotdir,
                             sample_locus_name,
                             "_",
                             "error_",
                             plotname,
                             "_",
                             mode,
                             ".png")
        } else if (type == 2) {
          filename <- paste0(plotdir,
                             paste0(gsub("[[:punct:]]", "", locus_id)),
                             "-",
                             sample_locus_name,
                             "_",
                             "error_",
                             plotname,
                             "_",
                             mode,
                             ".png")
        }

        plotmessage <- paste("Creating barplot No.", i)
        write_log(message = paste(plotmessage, "- filename:", filename),
                  logfilename = logfilename)

        dt <- data.table::data.table(
          "timepoint" = character(0),
          "value" = numeric(0),
          "regressiontype" = character(0)
        )

        # add relative error of corrected hyperbolic curve

        dt <- rbind(
          dt,
          cbind("timepoint" = "biased",
                "value" = round(
                  error_data[
                    get("Name") == vec_cal[i], get("relative_error_pre")
                  ],
                  3),
                "regressiontype" = "Raw")
        )
        dt <- rbind(
          dt,
          cbind("timepoint" = "corrected",
                "value" = round(
                  error_data[
                    get("Name") == vec_cal[i], get("relative_error")
                  ],
                  3),
                "regressiontype" = ifelse(
                  mode == "corrected_c",
                  "Corrected [Cubic]",
                  ifelse(
                    mode == "corrected_h",
                    "Corrected [Hyperbolic]",
                    "NA"))
          )
        )

        # set "Raw" as first level, to show corresponding bar
        # on the left of the plot
        dt[
          , ("regressiontype") := factor(get("regressiontype"),
                                         levels = c("Raw",
                                                    "Corrected [Cubic]",
                                                    "Corrected [Hyperbolic]")
          )
        ][
          , ("value") := round(as.numeric(as.character(get("value"))), 1)
        ]

        if ("Corrected [Cubic]" %in% dt[, get("regressiontype")]) {
          values <- c("#8491B4FF", "#E64B35FF")
        } else if ("Corrected [Hyperbolic]" %in% dt[, get("regressiontype")]) {
          values <- c("#8491B4FF", "#4DBBD5FF")
        }


        outplot <- ggplot2::ggplot(dt, ggplot2::aes_string(
          x = "regressiontype",
          y = "value",
          fill = "regressiontype")) +
          #% scale_fill_manual(
          #%   values = c("Cubic Regression" = "indianred1",
          #%              "Hyperbolic Regression" = "mediumspringgreen")) +
          ggplot2::geom_col() +
          ggplot2::geom_text(
            ggplot2::aes_string(
              label = "value",
              y = "value"),
            vjust = -1
          ) +
          ggplot2::ylab("% average relative error") +
          ggplot2::labs(
            title = paste0("Quantification Error: ", locus),
            subtitle = paste("CpG:", vec_cal[i]),
            fill = ggplot2::element_blank()
          ) +
          ggplot2::ylim(0, ylim_max) +
          ggplot2::scale_fill_manual(
            values = values
          ) +
          ggpubr::theme_pubr() +
          ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            text = ggplot2::element_text(size = plot_textsize)
          ) #,
        #% axis.ticks.x = element_blank(),
        #% axis.text.x = element_blank())
        #% print whole plot in return, otherwise it will fail

        ggplot2::ggsave(
          filename = filename,
          plot = outplot,
          device = "png",
          height = plot_height,
          width = plot_width,
          dpi = 600
        )
      })
    }, 1:length_vector)
  } else {
    write_log(
      message = "Error during creating bar plot; Names are not identical.",
      logfilename = logfilename)
  }
}

# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019 Lorenz Kapsner
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
                             headless = FALSE,
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

  Map(function(f) {
    plotname <- paste0(gsub("[[:punct:]]", "", rv$vec_cal[f]))

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

    # workaround to hide shiny-stuff, when going headless
    if (isFALSE(headless)) {
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = plotmessage, value = 0)

      # Increment the progress bar, and update the detail text.
      progress$inc(1 / 1, detail = paste("... working hard on plot",
                                         f,
                                         "of",
                                         length_vector))
    }

    # store plots to local temporary file
    create_plots(plotlist = plotlist_reg[[f]],
                 f = f,
                 rv = rv,
                 filename = filename,
                 logfilename = logfilename,
                 mode = mode,
                 minmax = minmax,
                 plot_height = plot_height,
                 plot_width = plot_width,
                 plot_textsize = plot_textsize)

  }, 1:length_vector)
}

create_plots <- function(plotlist,
                         f,
                         rv,
                         filename,
                         logfilename,
                         mode = NULL,
                         minmax,
                         plot_height = 5,
                         plot_width = 7.5,
                         plot_textsize = 16) {

  # get coefficients
  if (is.null(mode)) {
    # hyperbolic parameters
    coef_h <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]
    # cubic parameters
    coef_c <- rv$result_list[[rv$vec_cal[f]]][["Coef_cubic"]]
  } else if (mode == "corrected_h") {
    # hyperbolic parameters
    coef_h <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]
    # cubic parameters
    coef_c <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_cubic"]]
  } else if (mode == "corrected_c") {
    # hyperbolic parameters
    coef_h <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]
    # cubic parameters
    coef_c <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_cubic"]]
  }

  if (isFALSE(minmax)) {

    # create messages
    message <- paste0("# CpG-site: ", rv$vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", a =", coef_h$a,
                  ", b =", coef_h$b,
                  ", d =", coef_h$d)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", rv$vec_cal[f])
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
    message <- paste0("# CpG-site: ", rv$vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", y0 =", coef_h$y0,
                  ", y1 =", coef_h$y1,
                  ", m0 =", coef_h$m0,
                  ", m1 =", coef_h$m1)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", rv$vec_cal[f])
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
    width = plot_width
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
#' @export
#'
createbarerrorplots <- function(statstable_pre,
                                statstable_post,
                                rv,
                                type,
                                locus_id = NULL,
                                headless = FALSE,
                                plotdir,
                                logfilename,
                                mode = NULL,
                                plot_height = 5,
                                plot_width = 7.5,
                                plot_textsize = 16) {

  stats_pre <- statstable_pre[, c("Name", "relative_error"), with = F]
  stats_post <- statstable_post[, c("Name", "relative_error"), with = F]

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

    base::Map(function(i) {

      plotname <- paste0(gsub("[[:punct:]]", "", vec_cal[i]))

      if (type == 1) {
        filename <- paste0(plotdir,
                           rv$sample_locus_name,
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
                           rv$sample_locus_name,
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

      if (isFALSE(headless)) {
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive,
        # even if there's an error
        on.exit(progress$close())
        progress$set(message = plotmessage, value = 0)

        # Increment the progress bar, and update the detail text.
        progress$inc(
          1 / 1,
          detail = paste("... working hard on barplot",
                         i,
                         "of",
                         length_vector)
        )
      }

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
          , ("value") := as.numeric(as.character(get("value")))
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
        ggplot2::ylim(0, 100) +
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
        width = plot_width
      )

    }, 1:length_vector)
  } else {
    write_log(
      message = "Error during creating bar plot; Names are not identical.",
      logfilename = logfilename)
  }
}

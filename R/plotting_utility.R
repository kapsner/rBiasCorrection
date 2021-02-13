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
#' }
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

  future.apply::future_Map(
    f = function(f) {
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
    },
    1:length_vector,
    future.seed = TRUE
  )
}

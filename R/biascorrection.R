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


#' @title Correct PCR-Bias in Quantitative DNA Methylation Analyses.
#'
#' @description This function implements the algorithms described by
#'   Moskalev et. al in their article 'Correction of PCR-bias in quantitative
#'   DNA methylation studies by means of cubic polynomial regression',
#'   published 2011 in Nucleic acids research, Oxford University Press
#'   (\url{https://doi.org/10.1093/nar/gkr213}).
#'
#' @param experimental A character string. Path to the file containing the raw
#'   methylation values of the samples under investigation.
#' @param calibration A character string. In type 1 data (one locus in many
#'   samples, e.g. pyrosequencing data):
#'   Path to the file containing the raw methylation values of the calibration
#'   samples.
#'   In type 2 data (many loci in one sample, e.g. next-generationsequencing
#'   data or microarray data):
#'   Path to the folder that contains at least 4 calibration files (one file
#'   per calibration step). Please refere to the FAQ for more detailed
#'   information on the specific file requirements (\url{
#'   https://raw.githubusercontent.com/kapsner/PCRBiasCorrection/master/FAQ.md
#'   }).
#' @param samplelocusname A character string. In type 1 data: locus name -
#'   name of the gene locus under investigation. In type 2 data: sample name -
#'   name of the sample under investigation.
#' @param minmax A logical, indicating which equations are used for
#'   BiasCorrection (default: FALSE). If TRUE, equations are used that include
#'   the respective minima and maxima of the provided data.
#' @param correct_method A character string. Method used to correct the PCR-
#'   bias of the samples under investigation. One of "best" (default),
#'   "hyperbolic" or "cubic". If the method is set to "best" (short: "b"),
#'   the algorithm will automatically determine the best fitting type of
#'   regression for each CpG site based on \emph{selection_method} (by
#'   default: sum of squared errors, SSE,
#'   \url{https://en.wikipedia.org/wiki/Residual_sum_of_squares}). If the
#'   method is set to "hyperbolic" (short: "h") or "cubic" (short: "c"), the
#'   PCR-bias correction of all samples under investigation will be performed
#'   with the hyperbolic or the cubic regression respectively.
#' @param selection_method A character string. The method used to select the
#'   regression algorithm to correct the respective CpG site. This is by
#'   default the sum of squared errors ("SSE"). The second option is
#'   "RelError", which selects the regression method based on the theoretical
#'   relative error after correction. This metric is calculated by correcting
#'   the calibration data with both the hyperbolic regression and the cubic
#'   regression and using them again as input data to calculate the 'goodness
#'   of fit'-metrics.
#' @param type A single integer. Type of data to be corrected: either "1" (one
#'   locus in many samples, e.g. pyrosequencing data) or "2" (many loci in one
#'   sample, e.g. next-generation sequencing data or microarray data).
#' @param csvdir A character string. Directory to store the resulting tables.
#'   Default = "./csvdir". CAUTION: This directory will be newly created on
#'   every call of the function. Any preexisting files will be deleted without
#'   a warning.
#' @param plotdir A character string. Directory to store the resulting plots.
#'   Default = "./plotdir". CAUTION: This directory will be newly created on
#'   every call of the function. Any preexisting files will be deleted without
#'   a warning.
#' @param logfilename A character string. Path to a file to save the log
#'   messages. Default = "./log.txt"
#' @param seed A integer value. The seed used when solving the unknowns in the
#'   hyperbolic regression equation and the cubic regression equation.
#'   Important for reproducibility (default: 1234).
#' @param parallel A boolean. If TRUE (default), initializing
#'   `future::plan("multiprocess")` before running the code.
#'
#' @inheritParams createbarerrorplots
#'
#' @return This function is a wrapper around all of `rBiasCorrection`'s
#'   included functions. When executing it, it performs the whole workflow of
#'   bias correction and writes resulting csv-files and plots, as well as a
#'   log file to the local file system (the respective directories can be
#'   specified with the function arguments). The return-value is TRUE, if the
#'   correction of PCR measurement biases succeeds. If the correction fails,
#'   an error message is returned.
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' biascorrection("type1_experimentaldata.csv",
#'                "type1_calibrationdata.csv",
#'                samplelocusname = "BRAF")
#' }
#'
#' @export

biascorrection <- function(experimental,
                           calibration,
                           samplelocusname,
                           minmax = FALSE,
                           correct_method = "best",
                           selection_method = "SSE",
                           type = 1,
                           csvdir = "./csvdir",
                           plotdir = "./plotdir",
                           logfilename = "./log.txt",
                           plot_height = 5,
                           plot_width = 7.5,
                           plot_textsize = 16,
                           seed = 1234,
                           parallel = TRUE) {

  stopifnot(
    is.character(experimental),
    is.character(calibration),
    is.character(samplelocusname),
    is.logical(minmax),
    is.numeric(type),
    type == 1 || type == 2,
    is.character(correct_method),
    correct_method %in% c("best", "hyperbolic", "cubic", "b", "h", "c"),
    is.character(selection_method),
    selection_method %in% c("SSE", "RelError"),
    is.character(csvdir),
    is.character(plotdir),
    is.character(logfilename),
    is.numeric(plot_height),
    is.numeric(plot_width),
    is.numeric(plot_textsize),
    is.numeric(seed),
    is.logical(parallel)
  )

  # fix directories to work with all functions
  # therefore we need a "/" at the end of the dir-string
  plotdir <- gsub("([[:alnum:]])$", "\\1/", plotdir)
  csvdir <- gsub("([[:alnum:]])$", "\\1/", csvdir)

  # initialize some stuff
  on_start(plotdir, csvdir, logfilename, parallel)

  # initialize our list for reactive values
  rv <- list()

  # min-max hard coded to minmax = FALSE
  # (minmax = TRUE is experimental and under development)
  rv$minmax <- minmax
  rv$selection_method <- selection_method

  # save locusname
  rv$sample_locus_name <- handle_text_input(samplelocusname)

  # load data
  if (type == 1) {
    # experimental data
    rv$fileimport_experimental <- clean_dt(
      data.table::fread(experimental, header = T),
      "experimental",
      type = 1,
      logfilename = logfilename
    )[["dat"]]
    # write raw data to file
    write_csv(table = rv$fileimport_experimental,
              filename = paste0(csvdir, "raw_experimental_data.csv"))

    # calibration data
    cal_type_1 <- clean_dt(data.table::fread(calibration, header = T),
                           "calibration", type = 1,
                           logfilename = logfilename)
    rv$fileimport_calibration <- cal_type_1[["dat"]]
    # write raw data to file
    write_csv(table = rv$fileimport_calibration,
              filename = paste0(csvdir, "raw_calibration_data.csv"))

    # write names of columns to rv
    rv$vec_cal <- cal_type_1[["vec_cal"]]

    # calculate aggregated inputs
    rv$aggregated_experimental <- aggregated_input(
      datatable = rv$fileimport_experimental,
      description = "experimental",
      vec_cal = rv$vec_cal,
      type = 1
    )
    # write aggregated data to file
    write_csv(table = rv$aggregated_experimental,
              filename = paste0(csvdir, "aggregated_experimental_data.csv"))

    rv$aggregated_calibration <- aggregated_input(
      datatable = rv$fileimport_calibration,
      description = "calibration",
      vec_cal = rv$vec_cal
    )
    # write aggregated data to file
    write_csv(table = rv$aggregated_calibration,
              filename = paste0(csvdir, "aggregated_calibration_data.csv"))



  } else if (type == 2) {
    return(paste0("The correction of PCR measurement ",
                  "biases of this type of data is not implemented yet."))
  } else {
    return("ERROR. Please specify a valid type of data to correct (1 or 2).")
  }

  if (type == 1) {
    # calculate calibration curves
    # reconstruct parts from app_plottingUtility.R
    regression_results <- regression_utility(
      data = rv$fileimport_calibration,
      samplelocusname = rv$sample_locus_name,
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      headless = TRUE,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = seed
    )

    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

    # create calibration plots
    plotting_utility(
      data = rv$fileimport_calibration,
      plotlist_reg = plotlist_reg,
      type = 1,
      samplelocusname = rv$sample_locus_name,
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      headless = TRUE,
      plotdir = plotdir,
      logfilename = logfilename,
      minmax = rv$minmax,
      plot_height = plot_height,
      plot_width = plot_width,
      plot_textsize = plot_textsize
    )

    # save regression statistics to reactive value
    rv$reg_stats <- statistics_list(rv$result_list, minmax = rv$minmax)
    # write regression statistics to file
    write_csv(table = rv$reg_stats,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_regression_stats_",
                                gsub("\\-", "", substr(Sys.time(), 1, 10)),
                                "_",
                                gsub("\\:", "", substr(Sys.time(), 12, 16)),
                                ".csv")
    )

    # correct here calibration-data with hyperbolic and cubic regression
    # to get relative error of corrected data
    # hyperbolic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = F
                                    ][
                                      , ("better_model") := 0
                                      ]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_h <- solving_equations(datatable = rv$fileimport_calibration,
                                     regmethod = rv$choices_list,
                                     type = 1,
                                     rv = rv,
                                     mode = "corrected",
                                     logfilename = logfilename,
                                     minmax = rv$minmax)
    rv$fileimport_cal_corrected_h <- solved_eq_h[["results"]]
    colnames(rv$fileimport_cal_corrected_h) <- colnames(
      rv$fileimport_calibration
    )
    # write corrected calibration data to file
    write_csv(table = rv$fileimport_cal_corrected_h,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_corrected_calibrations_h_",
                                get_timestamp(),
                                ".csv"))

    # substitutions
    rv$substitutions_corrected_h <- solved_eq_h[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions_corrected_h) > 0) {
      write_csv(table = rv$substitutions_corrected_h,
                filename = paste0(csvdir,
                                  rv$sample_locus_name,
                                  "_substituted_corrected_h_",
                                  get_timestamp(),
                                  ".csv"))
    }

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      data = rv$fileimport_cal_corrected_h,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      headless = TRUE,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = seed
    )
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list_hyperbolic <- regression_results[["result_list"]]

    plotting_utility(
      data = rv$fileimport_cal_corrected_h,
      plotlist_reg = plotlist_reg,
      type = 1,
      samplelocusname = rv$sample_locus_name,
      locus_id = NULL,
      rv = rv,
      mode = "corrected_h",
      headless = TRUE,
      plotdir = plotdir,
      logfilename = logfilename,
      minmax = rv$minmax,
      plot_height = plot_height,
      plot_width = plot_width,
      plot_textsize = plot_textsize
    )

    # save regression statistics to reactive value
    rv$reg_stats_corrected_h <- statistics_list(
      resultlist = rv$result_list_hyperbolic,
      minmax = rv$minmax
    )
    # write regression statistics to file
    write_csv(table = rv$reg_stats_corrected_h,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_corrected_regression_stats_h_",
                                gsub("\\-", "", substr(Sys.time(), 1, 10)),
                                "_",
                                gsub("\\:", "", substr(Sys.time(), 12, 16)),
                                ".csv")
    )

    createbarerrorplots(
      statstable_pre = rv$reg_stats,
      statstable_post = rv$reg_stats_corrected_h,
      rv = rv,
      type = 1,
      locus_id = NULL,
      headless = TRUE,
      plotdir = plotdir,
      logfilename = logfilename,
      mode = "corrected_h",
      plot_height = plot_height,
      plot_width = plot_width,
      plot_textsize = plot_textsize
    )


    # cubic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = F
                                    ][
                                      , ("better_model") := 1
                                      ]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_c <- solving_equations(datatable = rv$fileimport_calibration,
                                     regmethod = rv$choices_list,
                                     type = 1,
                                     rv = rv,
                                     mode = "corrected",
                                     logfilename = logfilename,
                                     minmax = rv$minmax)
    rv$fileimport_cal_corrected_c <- solved_eq_c[["results"]]
    colnames(rv$fileimport_cal_corrected_c) <- colnames(
      rv$fileimport_calibration
    )
    # write corrected calibration data to file
    write_csv(table = rv$fileimport_cal_corrected_c,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_corrected_calibrations_c_",
                                get_timestamp(),
                                ".csv"))

    # substitutions
    rv$substitutions_corrected_c <- solved_eq_c[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions_corrected_c) > 0) {
      write_csv(table = rv$substitutions_corrected_c,
                filename = paste0(csvdir,
                                  rv$sample_locus_name,
                                  "_substituted_corrected_c_",
                                  get_timestamp(),
                                  ".csv")
      )
    }

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      data = rv$fileimport_cal_corrected_c,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      headless = TRUE,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = seed
    )
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list_cubic <- regression_results[["result_list"]]

    plotting_utility(
      data = rv$fileimport_experimental_corrected_c,
      plotlist_reg = plotlist_reg,
      type = 1,
      samplelocusname = rv$sample_locus_name,
      locus_id = NULL,
      rv = rv,
      mode = "corrected_c",
      headless = TRUE,
      plotdir = plotdir,
      logfilename = logfilename,
      minmax = rv$minmax,
      plot_height = plot_height,
      plot_width = plot_width,
      plot_textsize = plot_textsize
    )

    # save regression statistics to reactive value
    rv$reg_stats_corrected_c <- statistics_list(
      resultlist = rv$result_list_cubic,
      minmax = rv$minmax
    )
    # write regression statistics to file
    write_csv(table = rv$reg_stats_corrected_c,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_corrected_regression_stats_c_",
                                gsub("\\-", "", substr(Sys.time(), 1, 10)),
                                "_",
                                gsub("\\:", "", substr(Sys.time(), 12, 16)),
                                ".csv")
    )

    createbarerrorplots(
      statstable_pre = rv$reg_stats,
      statstable_post = rv$reg_stats_corrected_c,
      rv = rv,
      locus_id = NULL,
      type = 1,
      headless = TRUE,
      plotdir = plotdir,
      logfilename = logfilename,
      mode = "corrected_c",
      plot_height = plot_height,
      plot_width = plot_width,
      plot_textsize = plot_textsize
    )


    # now correct the real experimental data with the method chosen:
    # BiasCorrect experimental data with derived calibration curves
    if (correct_method %in% c("best", "b")) {
      # default selection of the model with the lower sse:
      rv$choices_list <- better_model(
        statstable_pre = rv$reg_stats,
        statstable_post_hyperbolic = rv$reg_stats_corrected_h,
        statstable_post_cubic = rv$reg_stats_corrected_c,
        selection_method = rv$selection_method
      )
    } else if (correct_method %in% c("hyperbolic", "h")) {
      rv$choices_list <- rv$reg_stats[, c("Name"), with = F
                                      ][
                                        , ("better_model") := 0
                                        ]
    } else if (correct_method %in% c("cubic", "c")) {
      rv$choices_list <- rv$reg_stats[, c("Name"), with = F
                                      ][
                                        , ("better_model") := 1
                                        ]
    }

    solved_eq <- solving_equations(
      datatable = rv$fileimport_experimental,
      regmethod = rv$choices_list,
      mode = NULL,
      type = 1,
      rv = rv,
      logfilename = logfilename,
      minmax = rv$minmax
    )
    rv$final_results <- solved_eq[["results"]]
    # write final results to csv
    write_csv(table = rv$final_results,
              filename = paste0(csvdir,
                                rv$sample_locus_name,
                                "_corrected_values_",
                                get_timestamp(),
                                ".csv"))
    rv$substitutions <- solved_eq[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions) > 0) {
      write_csv(table = rv$substitutions,
                filename = paste0(csvdir,
                                  rv$sample_locus_name,
                                  "_substituted_values_",
                                  get_timestamp(),
                                  ".csv"))
    }

  }

  return(TRUE)
}



# document datasets
#' @title example.data_experimental
#'
#' @description A list containing the experimental data ($dat) and the
#'   colnames ($vec_cal) needed by the algorithms to perform the bias
#'   correction.
#'
#' @name example.data_experimental
#' @docType data
#' @keywords data
NULL

#' @title example.data_calibration
#'
#' @description A list containing the calibration data ($dat) and the
#'   colnames ($vec_cal) needed by the algorithms to perform the bias
#'   correction.
#'
#' @name example.data_calibration
#' @docType data
#' @keywords data
NULL

#' @title example._plot.df_agg
#'
#' @description A data.table containing the aggregated
#'   calibration data for CpG site 1
#'   to create an example plot.
#'
#' @name example._plot.df_agg
#' @docType data
#' @keywords data
NULL

#' @title example._plot_coef_c
#'
#' @description A list containing exemplary coefficients
#'   of the cubic regression equation for CpG site 1
#'   to create an example plot.
#'
#' @name example._plot_coef_c
#' @docType data
#' @keywords data
NULL

#' @title example._plot_coef_h
#'
#' @description A list containing exemplary coefficients
#'   of the hyperbolic regression equation for CpG site 1
#'   to create an example plot.
#'
#' @name example._plot_coef_h
#' @docType data
#' @keywords data
NULL

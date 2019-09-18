# PCRBiasCorrection: Correct PCR-Bias in Quantitative DNA Methylation Analyses.
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


#' @title Correct PCR-Bias in Quantitative DNA Methylation Analyses.
#'
#' @description This function implements the algorithms described by Moskalev et. al in their
#'  article 'Correction of PCR-bias in quantitative DNA methylation studies by means of cubic polynomial regression',
#'  published 2011 in Nucleic acids research, Oxford University Press (\url{https://doi.org/10.1093/nar/gkr213}).
#'
#' @param experimental A character string. Path to the file containing the raw methylation
#'   values of the samples under investigation.
#' @param calibration A character string. In type 1 data (one locus in many samples, e.g. pyrosequencing data):
#'   Path to the file containing the raw methylation values of the calibration samples.
#'   In type 2 data (many loci in one sample, e.g. next-generation sequencing data or microarray data):
#'   Path to the folder that contains at least 4 calibration files (one file per calibration step).
#'   Please refere to the FAQ for more detailed information on the specific file requirements
#'   (\url{https://raw.githubusercontent.com/kapsner/PCRBiasCorrection/master/FAQ.md}).
#' @param samplelocusname A character string. In type 1 data: locus name - name of the gene locus under investigation.
#'   In type 2 data: sample name - name of the sample under investigation.
#' @param minmax A logical, indicating which equations are used for BiasCorrection (default: FALSE). If TRUE, equations are used
#'   that include the respective minima and maxima of the provided data.
#' @param correct_method A character string. Method used to correct the PCR bias of the samples under investigation.
#'   One of "best" (default), "hyperbolic" or "cubic". If the method is set to "best" (short: "b"), the algorithm will
#'   automatically determine the best fitting type of regression for each CpG site based on \emph{selection_method}
#'   (by default: sum of squared errors, SSE, \url{https://en.wikipedia.org/wiki/Residual_sum_of_squares}). If the method is set to
#'   "hyperbolic" (short: "h") or "cubic" (short: "c"), the PCR-bias correction of all samples under investigation will be performed
#'   with the hyperbolic or the cubic regression respectively.
#' @param selection_method A character string. The method used to select the regression algorithm to correct the
#'   respective CpG site. This is by default the sum of squared errors ("SSE"). The second option is "RelError",
#'   which selects the regression method based on the theoretical relative error after correction. This metric is
#'   calculated by correcting the calibration data with both the hyperbolic regression and the cubic regression and
#'   using them again as input data to calculate the 'goodness of fit'-metrics.
#' @param type A single integer. Type of data to be corrected: either "1" (one locus in many samples, e.g. pyrosequencing data)
#'   or "2" (many loci in one sample, e.g. next-generation sequencing data or microarray data).
#' @param csvdir A character string. Directory to store the resulting tables. Default = "./csvdir".
#'   CAUTION: This directory will be newly created on every call of the function. Any preexisting files
#'   will be deleted without a warning.
#' @param plotdir A character string. Directory to store the resulting plots. Default = "./plotdir".
#'   CAUTION: This directory will be newly created on every call of the function. Any preexisting files
#'   will be deleted without a warning.
#' @param logfilename A character string. Path to a file to save the log messages. Default = "./log.txt"
#'
#' @return TRUE, if the correction of PCR measurment biases succeeds. If the correction fails, an error message
#'   is returned.
#'
#' @import data.table
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' BiasCorrection("type1_experimentaldata.csv", "type1_calibrationdata.csv", samplelocusname = "BRAF")
#' }
#'
#' @export

BiasCorrection <- function(experimental, calibration, samplelocusname, minmax = FALSE, correct_method = "best", selection_method = "SSE", type = 1, csvdir = "./csvdir", plotdir = "./plotdir", logfilename = "./log.txt"){

  stopifnot(is.character(experimental),
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
            is.character(logfilename))

  # fix directories to work with all functions
  # therefore we need a "/" at the end of the dir-string
  plotdir <- gsub("([[:alnum:]])$", "\\1/", plotdir)
  csvdir <- gsub("([[:alnum:]])$", "\\1/", csvdir)

  # initialize some stuff
  onStart_(plotdir, csvdir, logfilename)


  # initialize our list for reactive values
  rv <- list()

  # min-max hard coded to minmax = FALSE (minmax = TRUE is experimental and under development)
  rv$minmax <- minmax
  rv$selection_method <- selection_method

  # save locusname
  rv$sampleLocusName <- handleTextInput_(samplelocusname)

  # load data
  if (type == 1){
    # experimental data
    rv$fileimportExp <- cleanDT_(data.table::fread(experimental, header = T), "experimental", type = 1, logfilename = logfilename)[["dat"]]
    # write raw data to file
    writeCSV_(rv$fileimportExp, paste0(csvdir, "raw_experimental_data.csv"))

    # calibration data
    cal_type_1 <- cleanDT_(data.table::fread(calibration, header = T), "calibration", type = 1, logfilename = logfilename)
    rv$fileimportCal <- cal_type_1[["dat"]]
    # write raw data to file
    writeCSV_(rv$fileimportCal, paste0(csvdir, "raw_calibration_data.csv"))

    # write names of columns to rv
    rv$vec_cal <- cal_type_1[["vec_cal"]]

  } else if (type == 2){
    return("The correction of PCR measurement Biases of this type of data is not implemented yet.")
  } else {
    return("ERROR. Please specify a valid type of data to correct (1 or 2).")
  }

  if (type == 1){
    # calculate calibration curves
    # reconstruct parts from app_plottingUtility.R
    regression_results <- regressionUtility_(rv$fileimportCal, rv$sampleLocusName, locus_id = NULL, rv = rv, mode = NULL, headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
    plotlistR <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

    # create calibration plots
    plottingUtility_(rv$fileimportCal, plotlistR, 1, rv$sampleLocusName, locus_id = NULL, rv = rv, mode = NULL, headless = TRUE, plotdir = plotdir, logfilename = logfilename, minmax = rv$minmax)

    # save regression statistics to reactive value
    rv$regStats <- statisticsList_(rv$result_list, minmax = rv$minmax)
    # write regression statistics to file
    writeCSV_(rv$regStats, #[,-(which(colnames(rv$regStats)=="better_model")), with=F],
             paste0(csvdir, "BC_regression_stats_", gsub("\\-", "", substr(Sys.time(), 1, 10)), "_",
                    gsub("\\:", "", substr(Sys.time(), 12, 16)), ".csv"))

    # correct here calibration-data with hyperbolic and cubic regression
    # to get relative error of corrected data
    # hyperbolic correction
    rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=0]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_h <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
    rv$fileimportCal_corrected_h <- solved_eq_h[["results"]]
    colnames(rv$fileimportCal_corrected_h) <- colnames(rv$fileimportCal)
    # write corrected calibration data to file
    writeCSV_(rv$fileimportCal_corrected_h,
              paste0(csvdir, "BC_corrected_calibrations_h_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))

    # substitutions
    rv$substitutions_corrected_h <- solved_eq_h[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions_corrected_h)>0){
      writeCSV_(rv$substitutions_corrected_h,
                paste0(csvdir, "BC_substituted_values_corrected_h_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))
    }

    # calculate new calibration curves from corrected calibration data
    regression_results <- regressionUtility_(rv$fileimportCal_corrected_h, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
    plotlistR <- regression_results[["plot_list"]]
    rv$result_list_hyperbolic <- regression_results[["result_list"]]

    plottingUtility_(rv$fileimportCal_corrected_h, plotlistR, type=1, samplelocusname=rv$sampleLocusName, locus_id = NULL, rv=rv, mode="corrected_h", headless = TRUE, plotdir = plotdir, logfilename = logfilename, minmax = rv$minmax)

    # save regression statistics to reactive value
    rv$regStats_corrected_h <- statisticsList_(rv$result_list_hyperbolic, minmax = rv$minmax)
    # write regression statistics to file
    writeCSV_(rv$regStats_corrected_h, #[,-(which(colnames(rv$regStats_corrected_h)=="better_model")), with=F],
              paste0(csvdir, "BC_corrected_regression_stats_h_", gsub("\\-", "", substr(Sys.time(), 1, 10)), "_",
                     gsub("\\:", "", substr(Sys.time(), 12, 16)), ".csv"))


    # for (i in rv$choices_list[,get("Name")]){
    #   rv$regStats_corrected_h[get("Name")==i,("better_model"):=rv$choices_list[get("Name")==i,as.integer(as.character(get("better_model")))]]
    # }
    createBarErrorPlots_(rv$regStats, rv$regStats_corrected_h, rv, type=1, headless = TRUE, plotdir = plotdir, logfilename = logfilename, mode = "corrected_h")


    # cubic correction
    rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=1]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_c <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
    rv$fileimportCal_corrected_c <- solved_eq_c[["results"]]
    colnames(rv$fileimportCal_corrected_c) <- colnames(rv$fileimportCal)
    # write corrected calibration data to file
    writeCSV_(rv$fileimportCal_corrected_c,
              paste0(csvdir, "BC_corrected_calibrations_c_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))

    # substitutions
    rv$substitutions_corrected_c <- solved_eq_c[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions_corrected_c)>0){
      writeCSV_(rv$substitutions_corrected_c,
                paste0(csvdir, "BC_substituted_values_corrected_c_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))
    }

    # calculate new calibration curves from corrected calibration data
    regression_results <- regressionUtility_(data=rv$fileimportCal_corrected_c, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
    plotlistR <- regression_results[["plot_list"]]
    rv$result_list_cubic <- regression_results[["result_list"]]

    plottingUtility_(rv$fileimportCal_corrected_c, plotlistR, type=1, samplelocusname=rv$sampleLocusName, locus_id = NULL, rv=rv, mode="corrected_c", headless = TRUE, plotdir = plotdir, logfilename = logfilename, minmax = rv$minmax)

    # save regression statistics to reactive value
    rv$regStats_corrected_c <- statisticsList_(rv$result_list_cubic, minmax = rv$minmax)
    # write regression statistics to file
    writeCSV_(rv$regStats_corrected_c, #[,-(which(colnames(rv$regStats_corrected_c)=="better_model")), with=F],
              paste0(csvdir, "BC_corrected_regression_stats_c_", gsub("\\-", "", substr(Sys.time(), 1, 10)), "_",
                     gsub("\\:", "", substr(Sys.time(), 12, 16)), ".csv"))


    # for (i in rv$choices_list[,get("Name")]){
    #   rv$regStats_corrected_c[get("Name")==i,("better_model"):=rv$choices_list[get("Name")==i,as.integer(as.character(get("better_model")))]]
    # }
    createBarErrorPlots_(rv$regStats, rv$regStats_corrected_c, rv, type=1, headless = TRUE, plotdir = plotdir, logfilename = logfilename, mode = "corrected_c")


    # now correct the real experimental data with the method chosen:
    # BiasCorrect experimental data with derived calibration curves
    if (correct_method %in% c("best", "b")){
      # default selection of the model with the lower sse:
      rv$choices_list <- betterModel(statstable_pre = rv$regStats,
                                     statstable_post_hyperbolic = rv$regStats_corrected_h,
                                     statstable_post_cubic = rv$regStats_corrected_c,
                                     selection_method = rv$selection_method)
    } else if (correct_method %in% c("hyperbolic", "h")){
      rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=0]
    } else if (correct_method %in% c("cubic", "c")){
      rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=1]
    }

    solved_eq <- solvingEquations_(rv$fileimportExp, rv$choices_list, type = 1, rv = rv, logfilename = logfilename, minmax = rv$minmax)
    rv$finalResults <- solved_eq[["results"]]
    # write final results to csv
    writeCSV_(rv$finalResults,
             paste0(csvdir, "BC_corrected_values_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))
    rv$substitutions <- solved_eq[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions)>0){
      writeCSV_(rv$substitutions,
               paste0(csvdir, "BC_substituted_values_", rv$sampleLocusName, "_", getTimestamp_(), ".csv"))
    }

  }

  return(TRUE)
}



# document datasets
#' @title example.data_experimental
#'
#' @description A list containing the experimental data ($dat) and the colnames ($vec_cal) needed by the algorithms to perform the bias correction.
#'
#' @name example.data_experimental
#' @docType data
#' @keywords data
NULL

#' @title example.data_calibration
#'
#' @description A list containing the calibration data ($dat) and the colnames ($vec_cal) needed by the algorithms to perform the bias correction.
#'
#' @name example.data_calibration
#' @docType data
#' @keywords data
NULL

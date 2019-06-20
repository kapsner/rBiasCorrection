# Correct PCR-Bias in Quantitative DNA Methylation Analyses.
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
#'   (\url{https://raw.githubusercontent.com/kapsner/BiasCorrector/master/FAQ.md}).
#' @param samplelocusname A character string. In type 1 data: locus name - name of the gene locus under investigation.
#'   In type 2 data: sample name - name of the sample under investigation.
#' @param type A single integer. Type of data to be corrected: either "1" (one locus in many samples, e.g. pyrosequencing data)
#'   or "2" (many loci in one sample, e.g. next-generation sequencing data or microarray data).
#' @param csvdir A character string. Path to a folder to store the resulting tables. Default = "./csvdir".
#' @param plotdir A character string. Path to a folder to store the resulting plots. Default = "./plotdir".
#' @param logfilename A character string. Path to a file to save the log messages. Default = "./log.txt"
#'
#' @return TRUE, if the correction of PCR measurment biases succeeds. If the correction fails, an error message
#'   is returned.
#'
#' @import data.table ggplot2 polynom ggpubr magrittr shiny
#'
#' @example
#' \dontrun{
#' BiasCorrection("type1_experimentaldata.csv", "type1_calibrationdata.csv", samplelocusname = "BRAF")
#' }
#'
#' @export

BiasCorrection <- function(experimental, calibration, samplelocusname, type=1, csvdir="./csvdir", plotdir="./plotdir", logfilename="./log.txt"){
  # fix directories to work with all functions
  # therefore we need a "/" at the end of the dir-string
  plotdir <- gsub("([[:alnum:]])$", "\\1/", plotdir)
  csvdir <- gsub("([[:alnum:]])$", "\\1/", csvdir)
  logfilename <<- logfilename

  # initialize some stuff
  onStart(plotdir, csvdir, logfilename)


  # initialize our list for reactive values
  rv <- list()

  # save locusname
  rv$sampleLocusName <- samplelocusname

  # load data
  if (type == 1){
    # experimental data
    rv$fileimportExp <- cleanDT(fread(experimental, header = T), "experimental", 1)[["dat"]]
    # write raw data to file
    writeCSV(rv$fileimportExp, paste0(csvdir, "raw_experimental_data.csv"))

    # calibration data
    cal_type_1 <- cleanDT(fread(calibration, header = T), "calibration", 1)
    rv$fileimportCal <- cal_type_1[["dat"]]
    # write raw data to file
    writeCSV(rv$fileimportCal, paste0(csvdir, "raw_calibration_data.csv"))

    # write names of columns to rv
    rv$vec_cal <- cal_type_1[["vec_cal"]]

  } else if (type == 2){
    print("The correction of PCR measurement Biases of this type of data is not implemented yet.")
  } else {
    return("ERROR. Please specify a valid type of data to correct (1 or 2).")
  }

  if (type == 1){
    # calculate calibration curves
    # reconstruct parts from app_plottingUtility.R
    regression_results <- regressionUtility(rv$fileimportCal, rv$sampleLocusName, locus_id = NULL, rv = rv, mode = NULL, headless = TRUE)
    plotlistR <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

    # create calibration plots
    plottingUtility(rv$fileimportCal, plotlistR, 1, rv$sampleLocusName, locus_id = NULL, rv = rv, mode = NULL, headless = TRUE, plotdir = plotdir)

    # save regression statistics to reactive value
    rv$regStats <- statisticsList(rv$result_list)
    # write regression statistics to file
    writeCSV(rv$regStats[,-(which(colnames(rv$regStats)=="better_model")), with=F],
             paste0(csvdir, "BC_regression_stats_", gsub("\\-", "", substr(Sys.time(), 1, 10)), "_",
                    gsub("\\:", "", substr(Sys.time(), 12, 16)), ".csv"))

    # BiasCorrect experimental data with derived calibration curves
    # default selection of the model with the lower sse:
    rv$choices_list <- rv$regStats[,c("Name", "better_model"),with=F]
    solved_eq <- solving_equations(rv$fileimportExp, rv$choices_list, type = 1, rv = rv)
    rv$finalResults <- solved_eq[["results"]]
    # write final results to csv
    writeCSV(rv$finalResults,
             paste0(csvdir, "BC_corrected_values_", rv$sampleLocusName, "_", getTimestamp(), ".csv"))
    rv$substitutions <- solved_eq[["substitutions"]]
    # write substitutions to csv (if existing)
    if (nrow(rv$substitutions)>0){
      writeCSV(rv$substitutions,
               paste0(csvdir, "BC_substituted_values_", rv$sampleLocusName, "_", getTimestamp(), ".csv"))
    }

    # extra stuff:
    # correct calibration data (to show corrected calibration curves)
    solved_eq2 <- solving_equations(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected")
    rv$fileimportCal_corrected <- solved_eq2[["results"]]
    colnames(rv$fileimportCal_corrected) <- colnames(rv$fileimportCal)
    # write corrected calibration data to file
    writeCSV(rv$fileimportCal_corrected,
             paste0(csvdir, "BC_corrected_calibrations_", rv$sampleLocusName, "_", getTimestamp(), ".csv"))

    # calculate new calibration curves from corrected calibration data
    regression_results <- regressionUtility(rv$fileimportCal_corrected, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE)
    plotlistR <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

    plottingUtility(rv$fileimportCal_corrected, plotlistR, type=1, samplelocusname=rv$sampleLocusName, locus_id = NULL, rv=rv, mode="corrected", headless = TRUE, plotdir = plotdir)

    # save regression statistics to reactive value
    rv$regStats_corrected <- statisticsList(rv$result_list)
    # write regression statistics to file
    writeCSV(rv$regStats_corrected[,-(which(colnames(rv$regStats_corrected)=="better_model")), with=F],
             paste0(csvdir, "BC_regression_stats_corrected_", gsub("\\-", "", substr(Sys.time(), 1, 10)), "_",
                    gsub("\\:", "", substr(Sys.time(), 12, 16)), ".csv"))


    for (i in rv$choices_list[,Name]){
      rv$regStats_corrected[Name==i,better_model:=rv$choices_list[Name==i,as.integer(as.character(better_model))]]
    }

    createBarErrorPlots(rv$regStats, rv$regStats_corrected, rv, type=1, headless = TRUE, plotdir = plotdir)
  }

  return(TRUE)
}

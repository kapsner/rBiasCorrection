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


#' @title type2FileReq helper function
#'
#' @description Function to check type 2 data file requirements.
#'
#' @export
#'
# check type 2 file requirements
type2FileReq_ <- function(filelist, rv){
  writeLog_("Entered 'type2FileReq'-Function")

  if (length(filelist) < 4){
    # error handling fileimport
    writeLog_("### ERROR ###\nLess than four calibration datasets were provided.\nAt least four distinct calibration steps are required to perform bias correction.")
    return("four")

    # go on, if there have been uploaded at least 4 files.
  } else {

    dims <- lapply(filelist, dim)

    # test if all files have the same dimension:
    if (length(unique(dims)) != 1){
      # error handling fileimport
      writeLog_("### ERROR ###\nThe dimensions of the provided calibration datasets are not equal.\nAll files have to have the same number of columns and rows.")
      return("dim")

      # if all files have the same dimension
    } else {

      # check here if all files have same column- and row-naming
      coln <- lapply(filelist, colnames)
      rown <- lapply(filelist, function(x){as.character(x[,locus_id])})

      if (length(unique(coln)) != 1 | length(unique(rown)) != 1){
        # error handling fileimport
        writeLog_("### ERROR ###\nAll files need to have the same rownames (locus ids) and columnnames (CpG sites).")
        return("naming")

      } else {

        # necessary for selectInput-choices
        #numberby <- 100/(length(rv$fileimportCal)-1)
        #choicesseq <- seq(from = 0, to = 100.00, by = numberby)
        pattern <- "(\\_CS\\d+(\\_\\d+)?(\\.csv|\\.CSV))$"
        calibr_steps <- data.table(name = character(), step = numeric())
        for (i in names(filelist)){
          message <- paste("Filename:", i)
          writeLog_(message)
          match <- regmatches(i, regexpr(pattern, i))
          if (length(match) < 1){
            writeLog_("### ERROR ###\nFilenaming of the calibration files must be done properly.\nEnd of filename must begin with '_CS' followd by a number, indicating the degree of true methylation.\nAs decimal seperator, '_' is required.")
            return("filename")
          } else {
            calibr_steps <- rbind(calibr_steps, data.table(name = i, step = as.numeric(gsub("\\_", ".", regmatches(match, regexpr("\\d+(\\_\\d+)?", match))))))
          }
        }
        calibr_steps <- calibr_steps[order(step, decreasing = F)]

        if (calibr_steps[,min(step)] < 0 | calibr_steps[,max(step)] > 100){
          writeLog_("### ERROR ###\nCalibration steps must be in range '0 <= calibration step <= 100'.")
          return("calibrange2")
        } else {

          # get unique gene names of first table (all tables must be equal, has been checked anywhere else??!)
          gene_names <- unique(filelist[[calibr_steps[1,name]]][,.(locus_id, CpG_count)])
          # get list of colnames
          col_names <- colnames(filelist[[calibr_steps[1,name]]])

          # check for duplicate locus_ids here
          if (sum(gene_names[,duplicated(locus_id)]) > 0){
            writeLog_("### ERROR ###\nPlease specify an equal number of CpG-sites for each gene locus.")
            return("inconsistency")
          } else {
            rv$calibr_steps <- calibr_steps
            return(TRUE)
          }
        }
      }
    }
  }
}

# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2022 Lorenz Kapsner
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


#' @title clean_dt helper function
#'
#' @description Internal function, that checks the formatting of imported
#'   files and prepares them to be applicable for the following pcr-bias
#'   correction steps.
#'
#' @param datatable A data.table object that contains either the experimental
#'   data or the calibration data.
#' @param description A character string, indicating if \code{datatable}
#'   contains either \emph{"calibration"} data or \emph{"experimental"} data.
#' @param type A single integer. Type of data to be corrected: either "1"
#'   (one locus in many samples, e.g. pyrosequencing data) or "2" (many loci
#'   in one sample, e.g. next-generation sequencing data or microarray data).
#' @param logfilename A character string. Path to the logfile to save the log
#'   messages.
#'
#' @return If a valid file is provided, the function returns a cleaned
#'   data.table, suited for BiasCorrection.
#'
#' @examples
#' logfilename <- paste0(tempdir(), "/log.txt")
#' cleaned_experimental <- clean_dt(
#'   datatable = rBiasCorrection::example.data_experimental$dat,
#'   description = "experimental",
#'   type = 1,
#'   logfilename = logfilename
#' )
#' dim(cleaned_experimental)
#' class(cleaned_experimental)
#'
#' @export
#'
clean_dt <- function(datatable, description, type, logfilename) {
  write_log(message = "Entered 'clean_dt'-Function",
            logfilename = logfilename)

  stopifnot(
    data.table::is.data.table(datatable),
    is.character(description),
    description %in% c("calibration", "experimental")
  )

  # workaround for vec_cal
  vec_cal <- NULL

  # remove all columns that only contain empty cells
  datatable <- Filter(function(x) !(all(x == "")), datatable)
  # remove all columns that only contain missings
  datatable <- Filter(function(x) !(all(is.na(x))), datatable)

  # load type 1 data
  if (type == "1") {
    message <- paste0("Importing data of type 1: One locus in many ",
                      "samples (e.g., pyrosequencing data)")
    write_log(message, logfilename)

    # rename cols, save colnames in vector
    if (description == "calibration") {
      message <- "got calibration data"
      write_log(message, logfilename)
      names(datatable)[1] <- "true_methylation"

    } else if (description == "experimental") {
      message <- "got experimental data"
      write_log(message, logfilename)
      names(datatable)[1] <- "sample_id"

    } else {
      message("### ERROR 18: wrong description ###")
      return(NULL)
    }

    # load type 2 data
  } else if (type == "2") {
    message <- paste0("Importing data of type 2: Many loci in one sample ",
                      "(e.g., next-gen seq or microarray data)")
    write_log(message, logfilename)

    # rename cols, save colnames in vector
    if (description == "calibration") {
      message <- "got calibration data"
      write_log(message, logfilename)
      names(datatable)[1] <- "locus_id"

    } else if (description == "experimental") {
      message <- "got experimental data"
      write_log(message, logfilename)
      names(datatable)[1] <- "locus_id"

    } else {
      message("### ERROR 36: wrong description ###")
      return(NULL)
    }
  } else {
    message("### ERROR 40: wrong type ###")
    return(NULL)
  }

  #% print(is.data.table(datatable))

  # all other columns are numeric
  vec <- names(datatable)
  datatable[, (vec) := lapply(.SD, function(x) {
    gsub(",", ".", x)
  }), .SDcols = vec]

  # fist column is factor
  if (description == "calibration" & type == "1") {
    # this is needed, because values here must be numeric
    result <- tryCatch({
      datatable[, (vec[1]) := factor(
        as.numeric(as.character(get(vec[1])))
      )]
    }, warning = function(w) {
      message("### ERROR 54: first column cannot be parsed to numeric ###")
      return(NULL)
    })

    if (is.null(result)) {
      return(NULL)
    }

  } else {
    datatable[, (vec[1]) := as.factor(get(vec[1]))]
  }

  # rest is numeric
  datatable[, (vec[-1]) := lapply(.SD, function(x) {
    as.numeric(as.character(x))
  }), .SDcols = vec[-1]]

  # sort datatable by first column and return it
  datatable <- datatable[order(datatable[[vec[1]]])]

  # remove empty rows
  datatable <- datatable[rowSums(datatable[, -1], na.rm = T) != 0, ]

  # some more dataprepration
  vec_cal <- names(datatable)[-1]

  # rowmeans are already in type 2 calibration data table
  # (from fileMerger-application!)
  datatable[, ("row_means") := rowMeans(
    datatable[, vec_cal, with = F], na.rm = T
    )]

  # make vec_cal global for type 1 data (many operations of
  # the app rely on vec_cal)
  if (type == "1") {
    vec_cal <- names(datatable)[-1]
  }

  # count number of CpGs in type 2 data
  if (type == "2") {
    datatable[, ("CpG_count") := rowSums(
      !is.na(datatable[, vec[-1], with = F])
      )]

    # requirements-check: does every repeated measurement of each locus id have
    # the same number of CpG-sites specified?
    if (sum(
      duplicated(
        unique(datatable[, get("CpG_count"), by = "locus_id"])$locus_id)
      ) > 0) {
      write_log(
        message = paste0("### ERROR ###\nThe data provided ",
                         "contains locus ids ",
                         "with heterogeneous counts of CpG-sites."),
        logfilename = logfilename
      )
      return(NULL)
    }

    if (description == "experimental") {
      # order experimental data by CpG-Count in decreasing order
      datatable <- datatable[order(get("CpG_count"), decreasing = T)]
    }
  }

  # check file requirements: missing values and remove rows
  # containing missing values
  if (type == "1") {
    #% if (isTRUE(any(is.na(datatable)))) {
    #%   before <- nrow(datatable)
    #%   datatable <- na.omit(datatable)
    #%   after <- nrow(datatable)
    #%   rv$omitnas <- before - after
    #%   message <- paste0("Deleted ",
    #%                     rv$omitnas,
    #%                     " row(s) containing missing values from '",
    #%                     description,
    #%                     " data'.")
    #%   write_log(message)
    #% }

    if (description == "calibration") {
      # type 1 data must have at least 4 calibration steps
      if (datatable[, nlevels(factor(get("true_methylation")))] < 4) {
        write_log(
          message = paste0("### ERROR ###\nThe data provided contains less ",
                           "than four calibration steps.\nAt least four ",
                           "distinct calibration steps are required to ",
                           "perform bias correction."),
          logfilename = logfilename
        )
        return(NULL)
      }
    }
  }
  return(list("dat" = datatable, "vec_cal" = vec_cal))
}

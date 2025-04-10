# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2025 Lorenz Kapsner
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

#' @title on_start helper function
#'
#' @description Internal function, that initializes plotdir,
#'   csvdir and logfilename.
#'
#' @param plotdir A character string. Path to the folder,
#'   where plots are saved.
#' @param csvdir A character string. Path to the folder,
#'   where resulting tables are saved.
#' @inheritParams clean_dt
#' @inheritParams biascorrection
#'
#' @return This function silently creates the directories`plotdir` and `csvdir`
#'   on the local filesystem and initializes the logfile, specified with
#'   `logfilename`. Furthermore, if `parallel = TRUE`, the `future`-backend is
#'   initialized.
#'
#' @examples
#' \donttest{ # runtime > 5 sec.
#' plotdir <- paste0(tempdir(), "/plots/")
#' csvdir <- paste0(tempdir(), "/csv/")
#' logfilename <- paste0(tempdir(), "/log.txt")
#' parallel <- FALSE
#'
#' on_start(plotdir, csvdir, logfilename, parallel)
#' }
#'
#' @seealso \link[future]{plan}
#'
#' @export
#'
on_start <- function(plotdir,
                     csvdir,
                     logfilename,
                     parallel) {

  if (dir.exists(plotdir)) {
    clean_up(plotdir, csvdir)
  }

  # create directories
  dir.create(plotdir)
  dir.create(csvdir)

  # initialize logfile here
  suppressMessages(suppressWarnings(file.create(logfilename)))

  if ("ggpubr" %in% utils::installed.packages()[, "Package"]) {
    options("ggpubr.exists" = TRUE)
  }

  if (isTRUE(parallel) && future::availableCores() > 1) {

    if (.Platform$OS.type == "unix") {
      write_log(
        message = "on_start: using future::plan(\"multicore\")",
        logfilename = logfilename
      )
      future::plan("multicore")
    } else {
      write_log(
        message = "on_start: using future::plan(\"multisession\")",
        logfilename = logfilename
      )
      future::plan("multisession")
    }
  } else {
    write_log(
      message = "on_start: using future::plan(\"sequential\")",
      logfilename = logfilename
    )
    future::plan("sequential")
  }
}

#' @title clean_up helper function
#'
#' @description Internal function to clean up directories.
#' @inheritParams on_start
#'
#' @return This function silently cleans up the current session and removes
#'   both, the `plotdir`- and the `csvdir`-folders. It furthermore resets
#'   the `future`-backend to plan = "sequential".
#'
#' @examples
#' plotdir <- paste0(tempdir(), "/plots/")
#' csvdir <- paste0(tempdir(), "/csv/")
#'
#' clean_up(plotdir, csvdir)
#'
#' @seealso \link[future]{plan}
#'
#' @export
#'
clean_up <- function(plotdir,
                     csvdir) {
  # activate sequential future
  suppressWarnings(future::plan("sequential"))

  # on session end, remove plots and and all other files from tempdir
  do.call(file.remove, list(list.files(plotdir, full.names = TRUE)))
  unlink(plotdir, recursive = TRUE)
  do.call(file.remove, list(list.files(csvdir, full.names = TRUE)))
  unlink(csvdir, recursive = TRUE)
}


#' @title write_log helper function
#'
#' @description Internal function to write log-messages to the file
#'   specified in logfilename.
#'
#' @param message A character string containing the log message.
#' @inheritParams clean_dt
#'
#' @return The function prints the logging message to the console and writes it
#'   to the local logfile, specified with `logfilename`.
#'
#' @examples
#' message <- "This is a logmessage"
#' logfilename <- paste0(tempdir(), "/log.txt")
#'
#' write_log(message, logfilename)
#'
#' @export
#'
# write log messages
write_log <- function(message, logfilename) {
  message(paste0("[", get_timestamp(), "]: ", message))
  message_out <- paste0("===========================================  \n",
                        "[Timestamp: ", get_timestamp(), "]  \n  \n",
                        message, "  \n  \n")
  write(message_out, file = logfilename, append = TRUE)
}


#' @title write_csv helper function
#'
#' @description Internal function to store the created tables in csv files.
#'
#' @param table A data.table object to store on the local file system
#' @param filename The file name (including the path) to store \code{table}.
#'
#' @return This function silently writes a `data.table` object to a CSV file.
#'
#' @examples
#' table <- data.table::data.table(
#'   a = stats::runif(1000),
#'   b = stats::runif(1000)
#' )
#'
#' write_csv(table, paste0(tempdir(), "/example.csv"))
#'
#' @seealso \link[data.table]{fwrite}
#'
#' @export
#'
# write csv files
write_csv <- function(table, filename) {
  return(data.table::fwrite(x = table,
                            file = filename,
                            row.names = FALSE,
                            sep = ",",
                            dec = ".",
                            eol = "\n"))
}

#' @title get_timestamp helper function
#'
#' @description Internal function to get the current timestamp to
#'   write it to filenames.
#'
#' @return This function takes no argument and returns a formatted timestamp
#'   of the current system time, which can be integrated e.g. into a filename.
#'
#' @examples
#' get_timestamp()
#'
#' @seealso \link{Sys.time}
#'
#' @export
# get_timestamp
get_timestamp <- function() {
  return(
    paste(gsub("\\-", "", substr(Sys.time(), 1, 10)),
          gsub("\\:", "", substr(Sys.time(), 12, 20)),
          sep = "_"))
}

# R-squared function
rsq <- function(true, fitted) {
  # https://en.wikipedia.org/wiki/Coefficient_of_determination
  return(stats::cor(true, fitted) ^ 2)
}


sdm <- function(vector) {
  return((vector - mean(vector))^2)
}

sse_tss <- function(datatable, fitted_values) {
  dat <- data.table::copy(datatable)

  # fitted values
  dat[, ("fitted") := fitted_values]

  # sum of squares between fitted and measured values
  dat[, ("CpG_fitted_diff") := get("CpG") - get("fitted")]
  dat[, ("squared_error") := get("CpG_fitted_diff")^2]

  # sum of squared errors = residual sum of squares
  sse <- as.numeric(dat[, sum(get("squared_error"), na.rm = TRUE)])

  # squared dist to mean
  dat[, ("squared_dist_mean") := sdm(get("fitted"))]

  # total sum of squares
  tss <- as.numeric(dat[, sum(get("squared_dist_mean"), na.rm = TRUE)])

  return(list(
    "sse" = sse,
    "tss" = tss
  ))
}


#' @title substitutions_create helper function
#'
#' @description Internal function to initialize a data.table object
#'   to store the substitutions.
#'
#' @return This function takes no argument and initializes an empty
#'   `data.table` to hold the substituted values, if substitutions occur during
#'   BiasCorrection.
#'
#' @examples
#' substitutions <- substitutions_create()
#' class(substitutions)
#'
#' @export
#'
# create substitutions dataframe
substitutions_create <- function() {
  substitutions <- data.table::data.table("id" = character(),
                                          "CpG_site" = character(),
                                          "corrected" = character(),
                                          "replacement" = character(),
                                          "regression" = character())
  return(substitutions)
}

#' @title handle_text_input helper function
#'
#' @description Internal function to remove punctuation and unneeded stuff
#'   from user inputs with regular expressions.
#'
#' @param textinput A character string with the textinput to perform these
#'   predefined regular expressions on.
#'
#' @return This function returns a cleaned up character string, limited to a
#'   maximum of 15 chars.
#'
#' @examples
#' textinput <- "This is a dirty! text."
#' handle_text_input(textinput)
#'
#' @export
#'
# handle user text inputs
handle_text_input <- function(textinput) {
  textinput <- gsub("[^[:alnum:]]", "", textinput)

  # max 15 chars:
  if (nchar(textinput) > 15) {
    textinput <- substr(textinput, 1, 15)
  }

  return(ifelse(nchar(textinput) > 0, textinput, "default"))
}


round_to_fifty <- function(max_err) {
  return(ceiling(max_err / 50) * 50)
}

testhelper_round_results_list <- function(x, dgts = 3) {
  list(
    Var = x$Var,
    relative_error = x$relative_error,
    SSE_hyper = x$SSE_hyper,
    Coef_hyper = lapply(x$Coef_hyper, round, digits = dgts),
    SSE_cubic = x$SSE_cubic,
    Coef_cubic = lapply(x$Coef_cubic, round, digits = dgts)
  )
}

testhelper_apply_robust_results_list <- function(res_list, dgts = 3) {
  sapply(
    X = res_list,
    FUN = testhelper_round_results_list,
    dgts = dgts,
    simplify = FALSE,
    USE.NAMES = TRUE
  )
}


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
#'
#' @export
#'
on_start <- function(plotdir,
                     csvdir,
                     logfilename) {

  if (dir.exists(plotdir)) {
    clean_up(plotdir, csvdir)
  }

  # create directories
  dir.create(plotdir)
  dir.create(csvdir)

  # initialize logfile here
  suppressMessages(suppressWarnings(file.create(logfilename)))
}

#' @title clean_up helper function
#'
#' @description Internal function to clean up directories.
#' @inheritParams on_start
#'
#' @export
#'
clean_up <- function(plotdir,
                     csvdir) {
  # on session end, remove plots and and all other files from tempdir
  do.call(file.remove, list(list.files(plotdir, full.names = TRUE)))
  unlink(plotdir, recursive = T)
  do.call(file.remove, list(list.files(csvdir, full.names = TRUE)))
  unlink(csvdir, recursive = T)
}


#' @title write_log helper function
#'
#' @description Internal function to write log-messages to the file
#'   specified in logfilename.
#'
#' @param message A character string containing the log message.
#' @inheritParams clean_dt
#'
#' @export
#'
# write log messages
write_log <- function(message, logfilename) {
  print(paste0("[", get_timestamp(), "]: ", message))
  message_out <- paste0("===========================================  \n",
                        "[Timestamp: ", get_timestamp(), "]  \n  \n",
                        message, "  \n  \n")
  write(message_out, file = logfilename, append = T)
}


#' @title write_csv helper function
#'
#' @description Internal function to store the created tables in csv files.
#'
#' @param table A data.table object to store on the local file system
#' @param filename The file name (including the path) to store \code{table}.
#'
#' @export
#'
# write csv files
write_csv <- function(table, filename) {
  return(data.table::fwrite(x = table,
                            file = filename,
                            row.names = F,
                            sep = ",",
                            dec = ".",
                            eol = "\n"))
}

#' @title get_timestamp helper function
#'
#' @description Internal function to get the current timestamp to
#'   write it to filenames.
#'
#' @export
#'
# get timestamp
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
  I((vector - mean(vector))^2)
}


#' @title substitutions_create helper function
#'
#' @description Internal function to initialize a data.table object
#'   to store the substitutions.
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

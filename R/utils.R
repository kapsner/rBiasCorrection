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

#' @title onStart helper function
#'
#' @description Initialization of plotdir, csvdir and logfilename
#'
#' @export
#'
onStart_ <- function(plotdir, csvdir, logfilename){

  if (dir.exists(plotdir)){
    cleanUp_(plotdir, csvdir)
  }

  # create directories
  dir.create(plotdir)
  dir.create(csvdir)

  # initialize logfile here
  suppressMessages(suppressWarnings(file.create(logfilename)))
}

#' @title cleanUp helper function
#'
#' @description Cleans up directories
#'
#' @export
#'
cleanUp_ <- function(plotdir, csvdir){
  # on session end, remove plots and and all other files from tempdir
  do.call(file.remove, list(list.files(plotdir, full.names = TRUE)))
  unlink(plotdir, recursive = T)
  do.call(file.remove, list(list.files(csvdir, full.names = TRUE)))
  unlink(csvdir, recursive = T)
}


#' @title writeLog helper function
#'
#' @description Writes log-messages to the file specified in logfilename
#'
#' @export
#'
# write log messages
writeLog_ <- function(message, logfilename){
  print(paste0("[", getTimestamp_(), "]: ", message))
  message_out <- paste0("===========================================  \n",
                        "[Timestamp: ", getTimestamp_(), "]  \n  \n",
                        message, "  \n  \n")
  write(message_out, file = logfilename, append = T)
}


#' @title writeCSV helper function
#'
#' @description Writes created tables to csv files
#'
#' @export
#'
# write csv files
writeCSV_ <- function(table, filename){
  return(data.table::fwrite(x = table,
                            file = filename,
                            row.names = F,
                            sep = ",",
                            dec = ".",
                            eol = "\n"))
}

#' @title getTimestamp helper function
#'
#' @description Gets the current timestamp to write it to filenames
#'
#' @export
#'
# get timestamp
getTimestamp_ <- function(){
  return(paste(gsub("\\-", "", substr(Sys.time(), 1, 10)), gsub("\\:", "", substr(Sys.time(), 12, 20)), sep="_"))
}

#' @title R-squared helper function
#'
#' @description Caclulates the Coefficient of determinition (R-squared,
#'   \url{https://en.wikipedia.org/wiki/Coefficient_of_determination}).
#'
# R-squared function
rsq <- function(true, fitted){
  return(stats::cor(true, fitted) ^ 2)
}


#' @title Squared distance to mean
#'
#' @description Calculates the squared distance to the mean in preparation for the calculation of R-squared
#'
sdm <- function(vector){
  I((vector-mean(vector))^2)
}


#' @title statisticsList helper function
#'
#' @description Formats the results_list into a data.table
#'
#' @export
#'
statisticsList_ <- function(resultlist){
  dt_list <- data.table::data.table("Name" = names(resultlist),
                                    "relative_error" = NA,
                                    "SSE_hyperbolic" = NA,
                                    "R2_hyperbolic" = NA,
                                    "b" = NA,
                                    "y0" = NA,
                                    "y1" = NA,
                                    "###" = NA,
                                    "SSE_cubic" = NA,
                                    "R2_cubic" = NA,
                                    "ax3" = NA,
                                    "bx2" = NA,
                                    "cx" = NA,
                                    "d" = NA)

  dt_list[, ("Name") := names(resultlist)]

  vec <- names(dt_list)[-1]
  dt_list[,(vec) := lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols = vec]

  for (i in names(resultlist)){
    dt_list[get("Name") == i, ("relative_error") := resultlist[[i]][["relative_error"]]
            ][
              get("Name") == i, ("SSE_hyperbolic") := resultlist[[i]][["SSE_hyper"]]
              ][
                get("Name") == i, ("R2_hyperbolic") := resultlist[[i]][["Coef_hyper"]][["R2"]]
                ][
                  get("Name") == i, ("b") := resultlist[[i]][["Coef_hyper"]][["b"]]
                  ][
                    get("Name") == i, ("y0") := resultlist[[i]][["Coef_hyper"]][["y0"]]
                    ][
                      get("Name") == i, ("y1") := resultlist[[i]][["Coef_hyper"]][["y1"]]
                      ][
                        get("Name") == i, ("SSE_cubic") := resultlist[[i]][["SSE_cubic"]]
                        ][
                          get("Name") == i, ("R2_cubic") := resultlist[[i]][["Coef_cubic"]][["R2"]]
                          ][
                            get("Name") == i, ("ax3") := resultlist[[i]][["Coef_cubic"]][["ax3"]]
                            ][
                              get("Name") == i, ("bx2") := resultlist[[i]][["Coef_cubic"]][["bx2"]]
                              ][
                                get("Name") == i, ("cx") := resultlist[[i]][["Coef_cubic"]][["cx"]]
                                ][
                                  get("Name") == i, ("d") := resultlist[[i]][["Coef_cubic"]][["d"]]
                                  ]
  }
  # mark the better model: 1 = cubic, 0 = hyperbolic
  dt_list[,("better_model") := ifelse(get("SSE_cubic") <= get("SSE_hyperbolic"), 1, 0)]

  return(dt_list)
}


#' @title substitutionsCreate helper function
#'
#' @description Function to initialize a data.table object to store the substitutions.
#'
#' @export
#'
# create substitutions dataframe
substitutionsCreate_ <- function(){
  substitutions <- data.table::data.table("id" = character(),
                                          "CpG_site" = character(),
                                          "corrected" = character(),
                                          "replacement" = character())
  return(substitutions)
}

#' @title handleTextInput helper function
#'
#' @description Function to remove punctuation and unneeded stuff from user inputs.
#'
#' @export
#'
# handle user text inputs
handleTextInput_ <- function(textinput){
  textinput <- gsub("[^[:alnum:]]", "", textinput)

  # max 15 chars:
  if (nchar(textinput) > 15){
    textinput <- substr(textinput, 1, 15)
  }

  return(ifelse(nchar(textinput) > 0, textinput, "default"))
}

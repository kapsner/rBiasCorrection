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

#' @title solving_equations helper function
#'
#' @description Internal function to solve the hyperbolic and cubic regression.
#'
#' @param regmethod A data.table object, with 2 columns, containing the names
#'   of the samples to correct (columns 1) and a binary variable
#'   \emph{better_model} that indicates, if the data should be corrected with
#'   the hyperbolic regression parameters (better_model = 0) or with the
#'   cubic regression parameters (better_model = 1).
#' @inheritParams clean_dt
#' @inheritParams regression_utility
#'
#' @return This function solves the equations of the hyperbolic and the
#'   cubic regression and returns the respectively interpolated values of the
#'   provided `datatable`.
#'
#' @examples
#' # define list object to save all data
#' rv <- list()
#' rv$minmax <- TRUE
#' rv$selection_method <- "RelError"
#' rv$sample_locus_name <- "Test"
#' rv$seed <- 1234
#'
#' # define logfilename
#' logfilename <- "log.txt"
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
#' # extract regression results
#' rv$result_list <- regression_results$result_list
#'
#' # get regression statistics
#' rv$reg_stats <- statistics_list(
#'   rv$result_list,
#'   minmax = TRUE
#' )
#'
#' # select the better model based on the sum of squared errrors ("SSE")
#' rv$choices_list <- better_model(
#'   statstable_pre = rv$reg_stats,
#'   selection_method = "SSE"
#' )
#'
#' # correct calibration data (to show corrected calibration curves)
#' solved_eq_h <- solving_equations(datatable = rv$fileimport_calibration,
#'                                  regmethod = rv$choices_list,
#'                                  type = 1,
#'                                  rv = rv,
#'                                  mode = "corrected",
#'                                  logfilename = logfilename,
#'                                  minmax = rv$minmax)
#' rv$fileimport_cal_corrected_h <- solved_eq_h$results
#' colnames(rv$fileimport_cal_corrected_h) <- colnames(
#'   rv$fileimport_calibration
#' )
#'
#' # cleanup
#' file.remove(logfilename)
#'
#'
#' @export
#'
# perform fitting of regressions to experimental data
solving_equations <- function(datatable,
                              regmethod,
                              type,
                              rv,
                              mode = NULL,
                              logfilename,
                              minmax) {
  write_log(
    message = "Entered 'solving_equations'-Function",
    logfilename = logfilename
  )

  # create data table to save results of substitutions
  substitutions <- substitutions_create()

  # get first colname
  first_colname <- colnames(datatable)[1]
  # create results dataframe and populate ids
  results <- data.table::data.table(
    "id" = datatable[, unique(get(first_colname))]
  )

  # loop through colnames aka. CpG-sites
  for (i in colnames(datatable)[-1]) {

    # initialize ouput-vector
    vector <- character()

    if (is.null(mode)) {
      df_agg_ex <- create_agg_df_exp(datatable, i, type)
    } else if (mode == "corrected") {
      df_agg_ex <- create_agg_df(datatable, i)
    }

    # if cubic regression has better sse-score (default), or
    # if user selects cubic regression for calculation manually in GUI
    if (regmethod[get("Name") == i, get("better_model")] == 1) {
      message <- paste("Solving cubic regression for", i)
      write_log(message = message, logfilename = logfilename)

      # loop through rows by samplenames
      for (j in as.vector(df_agg_ex[, get(first_colname)])) {
        msg1 <- paste("Samplename:", j)

        if (isFALSE(minmax)) {
          # get parameters
          ax3 <- rv$result_list[[i]][["Coef_cubic"]][["a"]]
          bx2 <- rv$result_list[[i]][["Coef_cubic"]][["b"]]
          cx <- rv$result_list[[i]][["Coef_cubic"]][["c"]]
          d <- rv$result_list[[i]][["Coef_cubic"]][["d"]]

          # this is the required form of the coefficients for
          # polynomial-function
          coe <- c(d - df_agg_ex[get(first_colname) == j, get("CpG")],
                   cx,
                   bx2,
                   ax3)

          #% print(paste(x_vec, class(x_vec)))
        } else if (isTRUE(minmax)) {
          a <- rv$result_list[[i]][["Coef_cubic"]][["a"]]
          b <- rv$result_list[[i]][["Coef_cubic"]][["b"]]
          y0 <- rv$result_list[[i]][["Coef_cubic"]][["y0"]]
          y1 <- rv$result_list[[i]][["Coef_cubic"]][["y1"]]
          m0 <- rv$result_list[[i]][["Coef_cubic"]][["m0"]]
          m1 <- rv$result_list[[i]][["Coef_cubic"]][["m1"]]

          # this is the required form of the coefficients for
          # polynomial-function
          c <- (((y1 - y0) / (m1 - m0)) - a * (m1 - m0)^2 - b * (m1 - m0))
          coe <- c(y0 - df_agg_ex[get(first_colname) == j, get("CpG")],
                   c,
                   b,
                   a)
        }

        print(coe)

        x_vec <- solve(polynom::polynomial(coe))               # polynom
        #% x_vec <- cubic(rev(coe))                    # RConics

        find_x <- list()

        for (m in x_vec) {
          #% print(m)
          if (grepl("(+|-)0i", as.character(m))) {   # muss am ende stehen!
            compl <- as.character(m)
            compl_out <- as.numeric((substr(compl, 1, nchar(compl) - 3)))
            find_x <- c(find_x, compl_out)
          } else {
            find_x <- c(find_x, m)
          }
        }

        if (isTRUE(minmax)) {
          # if minmax = TRUE and cubic = TRUE add m0, since we
          #% solved (x - m0) ^ j
          find_x <- lapply(find_x, function(x) {
            x + m0
          })
        }

        # generate checkpoint, to look, if fitting value has been found
        checkpoint <- FALSE
        # non-fitting numeric values
        nonfitting <- numeric()

        for (k in seq_len(length.out = length(find_x))) {
          if (class(find_x[[k]]) == "numeric") {
            if (find_x[[k]] >= 0 & find_x[[k]] <= 100) {
              # valid values must be 0 <= x <= 100
              msg2 <- "Root in between the borders! Added to results."
              vector <- c(vector, find_x[[k]])

              # break here, when first fitting value is found
              checkpoint <- TRUE
              write_log(
                message = paste0(msg1,
                                 "  \nRoot: ",
                                 round(find_x[[k]], 3),
                                 "  \n--> ",
                                 msg2),
                logfilename = logfilename
              )
              break

            } else {
              nonfitting <- c(nonfitting, find_x[[k]])
            }
          }
        }

        if (checkpoint == FALSE) {
          msg2 <- "## WARNING ##\nNo fitting root within the borders found."

          # if there no fitting numeric roots have been found, look,
          # if there are negative numeric roots
          # if there are nonfitting roots in range -10 - 0
          # (what we accept as "0")
          if (sum(nonfitting < 0 & nonfitting > -10) > 0) {
            nf <- nonfitting[nonfitting < 0 & nonfitting > -10]
            msg3 <- paste0("Negative numeric root found:  \nRoot: ",
                           round(nf, 3),
                           "  \n--> '-10 < root < 0' --> substitute 0")
            # if there are negative numeric roots, substitute value with "0",
            # which makes more sense!
            vector <- c(vector, 0)

            # store substitutions
            original <- as.character(nf)
            replacement <- "0"

          } else if (sum(nonfitting > 100 & nonfitting < 110) > 0) {
            nf <- nonfitting[nonfitting > 100 & nonfitting < 110]
            msg3 <- paste0("Positive numeric root found:  \nRoot: ",
                           round(nf, 3),
                           "  \n--> '100 < root < 110' --> substitute 100")
            # if there are positiv non fitting numeric roots in range 100 - 110
            # substitute value with "100", which makes more sense!
            vector <- c(vector, 100)

            # store substitutions
            original <- as.character(nf)
            replacement <- "100"

          } else {
            #if (is.null(mode)) {
            msg3 <- paste0("No fitting numeric roots within the ",
                           "borders found: substitute NA")
            vector <- c(vector, NA)

            # store substitutions
            original <- as.character(
              paste(round(nonfitting, 3),
                    collapse = ", ")
            )
            replacement <- "NA"
            #% } else {
            #%
            #%   print("Nonfitting")
            #%   print(nonfitting)
            #%
            #%   if (sum(nonfitting >= 110) > 0) {
            #%     nf <- nonfitting[nonfitting >= 110]
            #%     msg3 <- paste0("No fitting numeric roots within ",
            #%                    "the borders ",
            #%                    "found: since we are in mode == corrected, ",
            #%                    "substitute 100")
            #%     vector <- c(vector, 100)
            #%     replacement = "100"
            #%   } else if (sum(nonfitting <= -10) > 0) {
            #%     nf <- nonfitting[nonfitting <= -10]
            #%     msg3 <- paste0("No fitting numeric roots within ",
            #%                    "the borders ",
            #%                    "found: since we are in mode == corrected, ",
            #%                    "substitute 0")
            #%     vector <- c(vector, 0)
            #%     replacement = "0"
            #%   }
            #%   # store substitutions
            #%   original = as.character(paste(round(nonfitting, 3),
            #%                                 collapse = ", "))
            #% }

          }

          # substitutions only when correcting experimental --> No
          #if (is.null(mode)) {
          substitutions <- rbind(
            substitutions,
            data.table::data.table(
              "id" = j,
              "CpG_site" = i,
              "corrected" = original,
              "replacement" = replacement,
              "regression" = "cubic"
            ))
          #}

          write_log(
            message = paste0(msg1,
                             "  \n  \n",
                             msg2,
                             "  \n",
                             msg3),
            logfilename = logfilename
          )
        }
      }

      results <- results[, (paste0(i, "_c")) := as.numeric(
        as.character(
          get("vector")
        )
      )]

    } else if (regmethod[get("Name") == i, get("better_model")] == 0) {
      message <- paste("Solving hyperbolic regression for", i)
      write_log(message = message, logfilename = logfilename)

      for (j in as.vector(df_agg_ex[, get(first_colname)])) {
        msg1 <- paste("Samplename:", j)

        if (isFALSE(minmax)) {
          a <- rv$result_list[[i]][["Coef_hyper"]][["a"]]
          b <- rv$result_list[[i]][["Coef_hyper"]][["b"]]
          d <- rv$result_list[[i]][["Coef_hyper"]][["d"]]

          h_solv <- as.numeric(
            as.character(
              hyperbolic_eq_solved(
                y = df_agg_ex[get(first_colname) == j, get("CpG")],
                a = a,
                b = b,
                d = d
              )
            )
          )

        } else if (isTRUE(minmax)) {
          b <- rv$result_list[[i]][["Coef_hyper"]][["b"]]
          y0 <- rv$result_list[[i]][["Coef_hyper"]][["y0"]]
          y1 <- rv$result_list[[i]][["Coef_hyper"]][["y1"]]
          m0 <- rv$result_list[[i]][["Coef_hyper"]][["m0"]]
          m1 <- rv$result_list[[i]][["Coef_hyper"]][["m1"]]

          h_solv <- as.numeric(
            as.character(
              hyperbolic_eq_solved_minmax(
                y = df_agg_ex[get(first_colname) == j, get("CpG")],
                b = b,
                y0 = y0,
                y1 = y1,
                m0 = m0,
                m1 = m1)))
        }
        print(h_solv)

        if (h_solv >= 0 & h_solv <= 100) {
          msg2 <- "Root in between the borders! Added to results."
          vector <- c(vector, h_solv)
          write_log(
            message = paste0(msg1,
                             "  \nRoot: ",
                             round(h_solv, 3),
                             "  \n--> ",
                             msg2),
            logfilename = logfilename
          )

        } else {
          msg2 <- "## WARNING ##\nNo fitting root within the borders found."

          if (h_solv < 0 & h_solv > -10) {
            msg3 <- paste0("Negative numeric root found:  \nRoot: ",
                           round(h_solv, 3),
                           "  \n--> '-10 < root < 0' --> substitute 0")
            vector <- c(vector, 0)

            # store substitutions
            original <- as.character(h_solv)
            replacement <- "0"

          } else if (h_solv > 100 & h_solv < 110) {
            msg3 <- paste0("Positive numeric root found:  \nRoot: ",
                           round(h_solv, 3),
                           "  \n--> '100 < root < 110' --> substitute 100")
            vector <- c(vector, 100)

            # store substitutions
            original <- as.character(h_solv)
            replacement <- "100"

          } else {
            #if (is.null(mode)) {
            msg3 <- paste0("No fitting numeric roots within the borders ",
                           "found: substitute NA")
            vector <- c(vector, NA)

            # store substitutions
            original <- as.character(h_solv)
            replacement <- "NA"
            #% } else {
            #%   if (h_solv >= 110) {
            #%     msg3 <- paste0("No fitting numeric roots within ",
            #%                    "the borders ",
            #%                    "found: since we are in mode == corrected, ",
            #%                    "substitute 100")
            #%     vector <- c(vector, 100)
            #%     replacement = "100"
            #%   } else if (h_solv <= -10) {
            #%     msg3 <- paste0("No fitting numeric roots within ",
            #%                    "the borders ",
            #%                    "found: since we are in mode == corrected, ",
            #%                    "substitute 0")
            #%     vector <- c(vector, 0)
            #%     replacement = "0"
            #%   }
            #%   # store substitutions
            #%   original = as.character(h_solv)
            #% }


          }

          # substitutions only when correcting experimental --> No
          #if (is.null(mode)) {
          substitutions <- rbind(
            substitutions,
            data.table::data.table(
              "id" = j,
              "CpG_site" = i,
              "corrected" = original,
              "replacement" = replacement,
              "regression" = "hyperbolic"
            )
          )
          #}

          write_log(
            message = paste0(msg1,
                             "  \n  \n",
                             msg2,
                             "  \n",
                             msg3),
            logfilename = logfilename
          )
        }
      }

      # replace negative values, which are mathematically correct,
      # with "0", which makes more sense
      #% vector <- as.numeric(as.character(vector))
      #% vector <- ifelse(vector < 0,
      #%                  ifelse(vector > -10,
      #%                         0,
      #%                         NA),
      #%                  ifelse(vector <= 100,
      #%                         vector,
      #%                         ifelse(vector < 110,
      #%                                100,
      #%                                NA
      #%                         )
      #%                  )
      #% )

      # append output-vector to results
      results <- results[, (paste0(i, "_h")) := get("vector")]
    }
  }
  results[, (colnames(results)[-1]) := lapply(.SD, function(x) {
    as.numeric(as.character(x))
  }), .SDcols = colnames(results)[-1]]
  return(
    list("results" = results,
         "substitutions" = substitutions)
  )
}

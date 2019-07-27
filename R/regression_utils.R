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


# create aggregated datatable for calibration data
create_agg_df <- function(datatable, index){
  df <- datatable[,c("true_methylation", index), with = F]
  colnames(df)[2] <- "CpG"
  df[,("true_methylation") := as.numeric(as.character(get("true_methylation")))]
  return(df[, mean(get("CpG"), na.rm = T), by = "true_methylation"][,("CpG") := get("V1")][,("V1") := NULL])
}

# create aggregated datatable for experimental data
create_agg_df_exp <- function(datatable, index, type){
  if (type==1){
    df <- datatable[,c("sample_id", index), with = F]
    colnames(df)[2] <- "CpG"
    df <- df[, mean(get("CpG"), na.rm = T), by = "sample_id"][,("CpG") := get("V1")][,("V1") := NULL]
  } else if (type==2){
    df <- datatable[,c("locus_id", index), with = F]
    colnames(df)[2] <- "CpG"
    df <- df[, mean(get("CpG"), na.rm = T), by = "locus_id"][,("CpG") := get("V1")][,("V1") := NULL]
  }
  return(df)
}

#' @title solvingEquations helper function
#'
#' @description Solving hyperbolic and cubic regression.
#'
#' @export
#'
# perform fitting of regressions to experimental data
solvingEquations_ <- function(datatable, regmethod, type, rv, mode=NULL, logfilename, minmax){
  writeLog_("Entered 'solving_equations'-Function", logfilename)

  # create data table to save results of substitutions
  substitutions <- substitutionsCreate_()

  # get first colname
  first_colname <- colnames(datatable)[1]
  # create results dataframe and populate ids
  results <- data.table::data.table("id" = datatable[,unique(get(first_colname))])

  # loop through colnames aka. CpG-sites
  for (i in colnames(datatable)[-1]){

    # initialize ouput-vector
    vector <- character()

    if (is.null(mode)){
      df_agg_ex <- create_agg_df_exp(datatable, i, type)
    } else if (mode == "corrected"){
      df_agg_ex <- create_agg_df(datatable, i)
    }

    # if cubic regression has better sse-score (default), or
    # if user selects cubic regression for calculation manually in GUI
    if (regmethod[get("Name")==i,get("better_model")] == 1){
      message <- paste("Solving cubic regression for", i)
      writeLog_(message, logfilename)

      # get parameters
      ax3 <- rv$result_list[[i]][["Coef_cubic"]][["ax3"]]
      bx2 <- rv$result_list[[i]][["Coef_cubic"]][["bx2"]]
      cx <- rv$result_list[[i]][["Coef_cubic"]][["cx"]]
      d <- rv$result_list[[i]][["Coef_cubic"]][["d"]]

      # loop through rows by samplenames
      for (j in as.vector(df_agg_ex[,get(first_colname)])){
        msg1 <- paste("Samplename:", j)

        # this is the required form of the coefficients for polynomial-function
        coe <- c(d-df_agg_ex[get(first_colname)==j,get("CpG")], cx, bx2, ax3)
        print(coe)

        x_vec <- solve(polynom::polynomial(coe))               # polynom
        #x_vec <- cubic(rev(coe))                    # RConics

        #print(paste(x_vec, class(x_vec)))

        find_x <- list()

        for (m in x_vec){
          #print(m)
          if (grepl("(+|-)0i", as.character(m))){   # muss am ende stehen!
            compl <- as.character(m)
            compl_out <- as.numeric((substr(compl, 1, nchar(compl)-3)))
            find_x <- c(find_x, compl_out)
          } else {
            find_x <- c(find_x, m)
          }
        }

        # generate checkpoint, to look, if fitting value has been found
        checkpoint <- FALSE
        # non-fitting numeric values
        nonfitting <- numeric()

        for (k in 1:length(find_x)){
          if (class(find_x[[k]]) == "numeric"){
            if (find_x[[k]] >= 0 & find_x[[k]] <= 100){
              # valid values must be 0 <= x <= 100
              msg2 <- "Root in between the borders! Added to results."
              vector <- c(vector, find_x[[k]])

              # break here, when first fitting value is found
              checkpoint <- TRUE
              writeLog_(paste0(msg1, "  \nRoot: ", round(find_x[[k]], 3), "  \n--> ", msg2), logfilename)
              break

            } else {
              nonfitting <- c(nonfitting, find_x[[k]])
            }
          }
        }

        if (checkpoint == FALSE) {
          msg2 <- "## WARNING ##\nNo fitting root within the borders found."

          # if there no fitting numeric roots have been found, look, if there are negative numeric roots
          # if there are nonfitting roots in range -10 - 0 (what we accept as "0")
          if (sum(nonfitting < 0 & nonfitting > -10) > 0){
            nf <- nonfitting[nonfitting < 0 & nonfitting > -10]
            msg3 <- paste0("Negative numeric root found:  \nRoot: ", round(nf, 3), "  \n--> '-10 < root < 0' --> substitute 0")
            # if there are negative numeric roots, substitute value with "0", which makes more sense!
            vector <- c(vector, 0)

            # store substitutions
            original = as.character(nf)
            replacement = "0"

          } else if (sum(nonfitting > 100 & nonfitting < 110) > 0){
            nf <- nonfitting[nonfitting > 100 & nonfitting < 110]
            msg3 <- paste0("Positive numeric root found:  \nRoot: ", round(nf, 3), "  \n--> '100 < root < 110' --> substitute 100")
            # if there are positiv non fitting numeric roots in range 100 - 110
            # substitute value with "100", which makes more sense!
            vector <- c(vector, 100)

            # store substitutions
            original = as.character(nf)
            replacement = "100"

          } else {
            #if (is.null(mode)){
            msg3 <- "No fitting numeric roots within the borders found: substitute NA"
            vector <- c(vector, NA)

            # store substitutions
            original = as.character(paste(round(nonfitting, 3), collapse = ", "))
            replacement = "NA"
            # } else {
            #
            #   print("Nonfitting")
            #   print(nonfitting)
            #
            #   if (sum(nonfitting >= 110) > 0){
            #     nf <- nonfitting[nonfitting >= 110]
            #     msg3 <- "No fitting numeric roots within the borders found: since we are in mode==corrected,  substitute 100"
            #     vector <- c(vector, 100)
            #     replacement = "100"
            #   } else if (sum(nonfitting <= -10) > 0){
            #     nf <- nonfitting[nonfitting <= -10]
            #     msg3 <- "No fitting numeric roots within the borders found: since we are in mode==corrected,  substitute 0"
            #     vector <- c(vector, 0)
            #     replacement = "0"
            #   }
            #   # store substitutions
            #   original = as.character(paste(round(nonfitting, 3), collapse = ", "))
            # }

          }

          if (is.null(mode)){
            substitutions <- rbind(substitutions, data.table::data.table("id" = j,
                                                                         "CpG_site" = i,
                                                                         "corrected" = original,
                                                                         "replacement" = replacement))
          }

          writeLog_(paste0(msg1, "  \n  \n", msg2, "  \n", msg3), logfilename)
        }
      }

      results <- results[, (paste0(i, "_c")) := as.numeric(as.character(get("vector")))]

    } else if (regmethod[get("Name")==i,get("better_model")] == 0){
      message <- paste("Solving hyperbolic regression for", i)
      writeLog_(message, logfilename)

      if (isFALSE(minmax)){
        a <- rv$result_list[[i]][["Coef_hyper"]][["a"]]
        b <- rv$result_list[[i]][["Coef_hyper"]][["b"]]
        d <- rv$result_list[[i]][["Coef_hyper"]][["d"]]

      } else if (isTRUE(minmax)){
        b <- rv$result_list[[i]][["Coef_hyper"]][["b"]]
        y0 <- rv$result_list[[i]][["Coef_hyper"]][["y0"]]
        y1 <- rv$result_list[[i]][["Coef_hyper"]][["y1"]]
        m0 <- rv$result_list[[i]][["Coef_hyper"]][["m0"]]
        m1 <- rv$result_list[[i]][["Coef_hyper"]][["m1"]]
      }


      for (j in as.vector(df_agg_ex[,get(first_colname)])){
        msg1 <- paste("Samplename:", j)

        if (isFALSE(minmax)){
          h_solv <- as.numeric(as.character(hyperbolic_equation_solved(df_agg_ex[get(first_colname)==j,get("CpG")], a, b, d)))
        } else if (isTRUE(minmax)){
          h_solv <- as.numeric(as.character(hyperbolic_equation_solvedMinMax(df_agg_ex[get(first_colname)==j,get("CpG")], b, y0, y1, m0, m1)))
        }
        print(h_solv)

        if (h_solv >= 0 & h_solv <= 100){
          msg2 <- "Root in between the borders! Added to results."
          vector <- c(vector, h_solv)
          writeLog_(paste0(msg1, "  \nRoot: ", round(h_solv, 3), "  \n--> ", msg2), logfilename)

        } else {
          msg2 <- "## WARNING ##\nNo fitting root within the borders found."

          if (h_solv < 0 & h_solv > -10){
            msg3 <- paste0("Negative numeric root found:  \nRoot: ", round(h_solv, 3), "  \n--> '-10 < root < 0' --> substitute 0")
            vector <- c(vector, 0)

            # store substitutions
            original = as.character(h_solv)
            replacement = "0"

          } else if (h_solv > 100 & h_solv < 110){
            msg3 <- paste0("Positive numeric root found:  \nRoot: ", round(h_solv, 3), "  \n--> '100 < root < 110' --> substitute 100")
            vector <- c(vector, 100)

            # store substitutions
            original = as.character(h_solv)
            replacement = "100"

          } else {
            #if (is.null(mode)){
            msg3 <- "No fitting numeric roots within the borders found: substitute NA"
            vector <- c(vector, NA)

            # store substitutions
            original = as.character(h_solv)
            replacement = "NA"
            # } else {
            #   if (h_solv >= 110){
            #     msg3 <- "No fitting numeric roots within the borders found: since we are in mode==corrected,  substitute 100"
            #     vector <- c(vector, 100)
            #     replacement = "100"
            #   } else if (h_solv <= -10){
            #     msg3 <- "No fitting numeric roots within the borders found: since we are in mode==corrected,  substitute 0"
            #     vector <- c(vector, 0)
            #     replacement = "0"
            #   }
            #   # store substitutions
            #   original = as.character(h_solv)
            # }


          }

          if (is.null(mode)){
            substitutions <- rbind(substitutions, data.table::data.table("id" = j,
                                                                         "CpG_site" = i,
                                                                         "corrected" = original,
                                                                         "replacement" = replacement))
          }

          writeLog_(paste0(msg1, "  \n  \n", msg2, "  \n", msg3), logfilename)
        }
      }

      # replace negative values, which are mathematically correct, with "0", which makes more sense
      #vector <- as.numeric(as.character(vector))
      #vector <- ifelse(vector < 0, ifelse(vector > -10, 0, NA), ifelse(vector <= 100, vector, ifelse(vector < 110, 100, NA)))

      # append output-vector to results
      results <- results[, (paste0(i, "_h")) := get("vector")]
    }
  }
  results[,(colnames(results)[-1]):=lapply(.SD, function(x){as.numeric(as.character(x))}), .SDcols=colnames(results)[-1]]
  return(list("results" = results, "substitutions" = substitutions))
}

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


#' @title regressionUtility
#'
#' @description Internal function to carry out the regression calculations.
#'
#' @export
#'
regressionUtility_ <- function(data, samplelocusname, locus_id = NULL, rv, mode = NULL, headless = FALSE, logfilename, minmax = FALSE){

  if (!is.null(locus_id)){
    writeLog_(paste0("### Starting with regression calculations ###\n\nLocus ID: ", locus_id), logfilename)
  } else {
    writeLog_(paste0("### Starting with regression calculations ###"), logfilename)
  }


  # workaround to hide shiny-stuff, when going headless
  if (isFALSE(headless)){
    # for plotting: basic idea and some code snippets from:
    # https://gist.github.com/wch/5436415/
    regression <- shiny::reactive({
      regression_type1(data, rv$vec_cal, mode, logfilename)
    })

    shiny::withProgress(message = "Calculating calibration curves", value = 0, {
      shiny::incProgress(1/1, detail = "... working on calculations ...")
      # calculate results (if this is run here, j must be resetted)
      regression_results <- regression()
    })
  } else {
    regression_results <- regression_type1(data, rv$vec_cal, mode, logfilename, minmax = minmax)
  }

  return(list("plot_list" = regression_results[["plot_list"]],
              "result_list" = regression_results[["result_list"]]))
}

regression_type1 <- function(datatable, vec_cal, mode=NULL, logfilename, minmax){
  writeLog_("Entered 'regression_type1'-Function", logfilename)

  # result_list
  result_list <- list()

  plot.listR <- list()

  for (i in 1:length(vec_cal)){
    message <- paste0("# CpG-site: ", vec_cal[i])
    writeLog_(message, logfilename)
    df_agg <- stats::na.omit(create_agg_df(datatable, vec_cal[i]))

    print(df_agg)
    writeLog_(paste("Logging df_agg:", vec_cal[i]), logfilename)
    writeLog_(df_agg, logfilename)

    result_list[[vec_cal[i]]] <- hyperbolic_regression(df_agg, vec_cal[i], logfilename, minmax = minmax)
    # append result_list
    result_list[[vec_cal[i]]] <- c(result_list[[vec_cal[i]]], cubic_regression(df_agg, vec_cal[i], logfilename))

    if (is.null(mode)){
      custom_ylab <- "% apparent methylation after PCR"
    } else if (mode == "corrected"){
      custom_ylab <- "% methylation after BiasCorrection"
    }

    lb1 <- c(paste0(" Hyperbolic: R\u00B2=", round(result_list[[vec_cal[i]]]$Coef_hyper$R2, 2)),
             paste0(" Cubic: R\u00B2=", round(result_list[[vec_cal[i]]]$Coef_cubic$R2, 2))
    )


    gdat <- df_agg[
      ,("true_methylation"):=as.numeric(as.character(get("true_methylation")))
      ][
        ,("CpG"):=as.numeric(as.character(get("CpG")))
      ]

    p <- ggplot2::ggplot(data=gdat, ggplot2::aes_string(x = "true_methylation", y = "CpG")) +
      ggplot2::geom_point() +
      ggplot2::ylab(custom_ylab) +
      ggplot2::xlab("% actual methylation") +
      ggplot2::ggtitle(paste("CpG-site:", vec_cal[i])) +
      ggplot2:: geom_text(data = data.frame(),
                          ggplot2::aes(x=-Inf, y=c(max(gdat$CpG), 0.95*max(gdat$CpG)), hjust=0, vjust = 1),
                          label = lb1,
                          size = 3.5,
                          parse = F)
    plot.listR[[i]] <- p
  }
  return(list("plot_list" = plot.listR, "result_list" = result_list))
}

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

#' @title plottingUtility helper function
#'
#' @description Function to carry out the plotting of the calibrations curves.
#'
#' @export
#'
# plotting utility
plottingUtility_ <- function(data, plotlistR, type, samplelocusname, locus_id = NULL, rv, mode=NULL, headless=FALSE, plotdir, logfilename){

  if (!is.null(locus_id)){
    writeLog_(paste0("### Starting with plotting ###\n\nLocus ID: ", locus_id), logfilename)
  } else {
    writeLog_(paste0("### Starting with plotting ###"), logfilename)
  }

  # get number of CpG-sites
  length_vector <- length(rv$vec_cal)

  Map(function(f) {
    plotname <- paste0(gsub("[[:punct:]]", "", rv$vec_cal[f]))

    # filename-suffix
    fn_suffix <- ifelse(is.null(mode), "", paste0("_", mode))
    # message-suffix
    msg_suffix <- ifelse(is.null(mode), "", ifelse(mode == "corrected_h", "BiasCorrected (hyperbolic)", "BiasCorrected (cubic)"))

    # filname of temporary plot
    if (type == 1){
      filename <- paste0(plotdir, samplelocusname, "_", plotname, fn_suffix, ".png")
      plotmessage <- paste0("Creating ", msg_suffix, " plot No. ", f)
    } else if (type == 2){
      filename <- paste0(plotdir, locus_id, "-", samplelocusname, "_", plotname, fn_suffix, ".png")
      plotmessage <- paste0("Locus ID: ", locus_id, " --> Creating ", msg_suffix, " plot No. ", f)
    }

    writeLog_(paste(plotmessage, "- filename:", filename), logfilename)

    # workaround to hide shiny-stuff, when going headless
    if (isFALSE(headless)){
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = plotmessage, value = 0)

      # Increment the progress bar, and update the detail text.
      progress$inc(1/1, detail = paste("... working hard on plot", f, "of", length_vector))
    }

    # store plots to local temporary file
    createPlots(plotlistR[[f]], f, rv, filename, logfilename, mode)

  }, 1:length_vector)
}

createPlots <- function(plotlist, f, rv, filename, logfilename, mode = NULL){
  shiny::plotPNG({

    if (is.null(mode)){
      # hyperbolic parameters
      b <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]$b
      y0 <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]$y0
      y1 <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]$y1
      m0 <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]$m0
      m1 <- rv$result_list[[rv$vec_cal[f]]][["Coef_hyper"]]$m1

      # cubic parameters
      c <- sapply(rv$result_list[[rv$vec_cal[f]]][["Coef_cubic"]], unlist)[c(4:1)]
      #c <- sapply(rv$result_list[[rv$vec_cal[i]]][["Coef_cubic"]], `[`)[c(4:1)]
    } else if (mode == "corrected_h"){
      # hyperbolic parameters
      b <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]$b
      y0 <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]$y0
      y1 <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]$y1
      m0 <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]$m0
      m1 <- rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_hyper"]]$m1

      # cubic parameters
      c <- sapply(rv$result_list_hyperbolic[[rv$vec_cal[f]]][["Coef_cubic"]], unlist)[c(4:1)]
      #c <- sapply(rv$result_list[[rv$vec_cal[i]]][["Coef_cubic"]], `[`)[c(4:1)]
    } else if (mode == "corrected_c"){
      # hyperbolic parameters
      b <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]$b
      y0 <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]$y0
      y1 <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]$y1
      m0 <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]$m0
      m1 <- rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_hyper"]]$m1

      # cubic parameters
      c <- sapply(rv$result_list_cubic[[rv$vec_cal[f]]][["Coef_cubic"]], unlist)[c(4:1)]
      #c <- sapply(rv$result_list[[rv$vec_cal[i]]][["Coef_cubic"]], `[`)[c(4:1)]
    }



    # create messages
    message <- paste0("# CpG-site: ", rv$vec_cal[f])
    msg2 <- paste("Using bias_weight =", b, ", y0 =", y0, ", y1 =", y1)
    writeLog_(paste(message, msg2, sep="\n"), logfilename)

    message <- paste0("# CpG-site: ", rv$vec_cal[f])
    msg2 <- paste("Using c =", paste(c, collapse = ", "))
    writeLog_(paste0(message, msg2, sep="\n"), logfilename)

    return(print(plotlist +
                   ggplot2::stat_function(fun = hyperbolic_equation, args = list(b=b, y0=y0, y1=y1, m0=m0, m1=m1), geom = "line", ggplot2::aes(color = "Hyperbolic Regression"), size=1.06) +
                   ggplot2::stat_function(fun = cubic_equation, args = list(c=c), geom = "line", ggplot2::aes(color = "Cubic Regression"), size=1.06) +
                   ggplot2::geom_line(ggplot2::aes(x=plotlist$data$true_methylation, y=plotlist$data$true_methylation, color = "unbiased"), linetype="dashed", size=1.04) +
                   ggplot2::labs(color = ggplot2::element_blank()) +
                   ggplot2::scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF"),
                                               labels = c("Cubic Regression", "Hyperbolic Regression", "unbiased")) +
                   #scale_colour_manual("Regression:", values = c(Cubic = "indianred1", Hyperbolic = "mediumspringgreen", unbiased = "lightblue")) +
                   ggpubr::theme_pubr() +
                   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    ))
  },
  filename = filename,
  height = 400,
  width = 450)
}


#' @title createBarErrorPlots helper function
#'
#' @description Function to create relative-error bar plots.
#'
#' @export
#'
createBarErrorPlots_ <- function(statstable_pre, statstable_post, rv, type, b=NULL, headless = FALSE, plotdir, logfilename, mode = NULL){

  stats_pre <- statstable_pre[,c("Name", "relative_error", "better_model"),with=F]
  stats_post <- statstable_post[,c("Name", "relative_error", "better_model"),with=F]

  error_data <- merge(stats_post, stats_pre, by="Name", sort=F, suffixes=c("", "_pre"))

  # Test if names are eqal
  if (identical(stats_pre[,get("Name")], stats_post[,get("Name")])){

    vec_cal <- stats_pre[,get("Name")]
    length_vector <- length(vec_cal)

    base::Map(function(i) {

      plotname <- paste0(gsub("[[:punct:]]", "", vec_cal[i]))

      if (type == 1){
        filename <- paste0(plotdir, "Errorplot_", rv$sampleLocusName, "_", plotname, "_", mode, ".png")
      } else if (type == 2){
        filename <- paste0(plotdir, "Errorplot_", paste0(gsub("[[:punct:]]", "", b)), "-", rv$sampleLocusName, "_", plotname, "_", mode, ".png")
      }

      plotmessage <- paste("Creating barplot No.", i)
      writeLog_(paste(plotmessage, "- filename:", filename), logfilename)

      dt <- data.table::data.table("timepoint" = character(0), "value" = numeric(0), "regressiontype" = character(0))

      if (isFALSE(headless)){
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = plotmessage, value = 0)

        # Increment the progress bar, and update the detail text.
        progress$inc(1/1, detail = paste("... working hard on barplot", i, "of", length_vector))
      }

      # add relative error of corrected hyperbolic curve

      dt <- rbind(dt, cbind("timepoint"="biased", "value" = round(error_data[get("Name")==vec_cal[i],get("relative_error_pre")], 3),
                            "regressiontype" = "Uncorrected"))
      dt <- rbind(dt, cbind("timepoint"="corrected", "value" = round(error_data[get("Name")==vec_cal[i],get("relative_error")], 3),
                            "regressiontype" = ifelse(mode == "corrected_c", "Corrected [Cubic Regression]",
                                                      ifelse(mode == "corrected_h", "Corrected [Hyperbolic Regression]", "NA"))))

      # set "Raw" as first level, to show corresponding bar on the left of the plot
      dt[
        ,("regressiontype") := factor(get("regressiontype"), levels = c("Uncorrected", "Corrected [Cubic Regression]", "Corrected [Hyperbolic Regression]"))
        ][
          ,("value"):=as.numeric(as.character(get("value")))
        ]

      if ("Corrected [Cubic Regression]" %in% dt[,get("regressiontype")]){
        values <- c("#8491B4FF", "#E64B35FF")
      } else if ("Corrected [Hyperbolic Regression]" %in% dt[,get("regressiontype")]){
        values <- c("#8491B4FF", "#4DBBD5FF")
      }

      shiny::plotPNG({
        p <- ggplot2::ggplot(dt, ggplot2::aes_string(x = "regressiontype", y="value", fill="regressiontype")) +
          #scale_fill_manual(values = c("Cubic Regression" = "indianred1", "Hyperbolic Regression" = "mediumspringgreen")) +
          ggplot2::geom_col()+
          ggplot2::geom_text(ggplot2::aes_string(label = "value", y="value"),  vjust = 3) +
          ggplot2::ylab("% average relative error") +
          ggplot2::labs(title = paste0("Quantification Error ", vec_cal[i]), fill = ggplot2::element_blank()) +
          ggplot2::scale_fill_manual(values = values) +
          ggpubr::theme_pubr() +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         legend.position = "none",
                         plot.title = ggplot2::element_text(hjust = 0.5)) #,
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank())
        # print whole plot in return, otherwise it will fail
        return(print(p))
      },
      filename = filename,
      height = 400,
      width = 450)
    }, 1:length_vector)
  } else {
    writeLog_("Error during creating bar plot; Names are not identical.", logfilename)
  }
}

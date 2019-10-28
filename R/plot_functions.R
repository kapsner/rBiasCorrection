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

#' @title calibration_plot helper function
#'
#' @description Internal function to carry out the plotting of the
#'   calibrations curves.
#'
#' @param plotlist A ggplot object contating a regression plot without
#'   regression curves (output of \code{regression_utility()}).
#' @param coef_hyper A list containing the regression parameters of
#'   the hyperbolic regression equation.
#' @param coef_cubic A list containing the regression parameters of
#'   the cubic regression equation.
#' @inheritParams regression_utility
#' @inheritParams createbarerrorplots
#'
#' @export
#'
# calibration_plot
calibration_plot <- function(plotlist,
                             coef_hyper,
                             coef_cubic,
                             plot_textsize,
                             minmax) {

  if (isFALSE(minmax)) {
    p <- plotlist +
      ggplot2::stat_function(
        fun = hyperbolic_eq,
        args = list(a = coef_hyper$a,
                    b = coef_hyper$b,
                    d = coef_hyper$d),
        geom = "line",
        ggplot2::aes(
          color = "Hyperbolic Regression"),
        size = 1.08
      ) +
      ggplot2::stat_function(
        fun = cubic_eq,
        args = list(a = coef_cubic$a,
                    b = coef_cubic$b,
                    c = coef_cubic$c,
                    d = coef_cubic$d),
        geom = "line",
        ggplot2::aes(color = "Cubic Regression"),
        size = 1.08
      )

  } else if (isTRUE(minmax)) {

    p <- plotlist +
      ggplot2::stat_function(
        fun = hyperbolic_eq_minmax,
        args = list(b = coef_hyper$b,
                    y0 = coef_hyper$y0,
                    y1 = coef_hyper$y1,
                    m0 = coef_hyper$m0,
                    m1 = coef_hyper$m1),
        geom = "line",
        ggplot2::aes(
          color = "Hyperbolic Regression"),
        size = 1.08
      ) +
      ggplot2::stat_function(
        fun = cubic_eq_minmax,
        args = list(a = coef_cubic$a,
                    b = coef_cubic$b,
                    y0 = coef_cubic$y0,
                    y1 = coef_cubic$y1,
                    m0 = coef_cubic$m0,
                    m1 = coef_cubic$m1),
        geom = "line", ggplot2::aes(
          color = "Cubic Regression"),
        size = 1.08
      )
  }

  outplot <- p +
    ggplot2::geom_line(
      ggplot2::aes(
        x = plotlist$data$true_methylation,
        y = plotlist$data$true_methylation,
        color = "unbiased"),
      linetype = "dashed",
      size = 1.08
    ) +
    ggplot2::labs(
      color = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(
      values = c("#E64B35FF",
                 "#4DBBD5FF",
                 "#00A087FF"),
      labels = c("Cubic Regression",
                 "Hyperbolic Regression",
                 "unbiased")
    ) +
    #% scale_colour_manual("Regression:",
    #%                     values = c(Cubic = "indianred1",
    #%                                Hyperbolic = "mediumspringgreen",
    #%                                unbiased = "lightblue")) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(size = plot_textsize)
    )

  return(outplot)
}


#' @title create_exampleplot helper function
#'
#' @description Internal function to create an example plot.
#'
#' @param data A data.table containing the aggregated
#'   calibration data.
#' @param filename A character. The filename, where to
#'   store the resulting example plot.
#'
#' @inheritParams calibration_plot
#' @inheritParams regression_utility
#' @inheritParams createbarerrorplots
#'
#' @export
#'
# create_exampleplot
create_exampleplot <- function(data,
                               coef_hyper,
                               coef_cubic,
                               plot_height,
                               plot_width,
                               plot_textsize,
                               filename) {

  lb1 <- c(paste0("  R\u00B2: \n  Hyperbolic = ",
                  round(coef_hyper$R2, 2),
                  "\n  Cubic = ",
                  round(coef_cubic$R2, 2)), "")

  # create base plot
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(
                         x = "true_methylation",
                         y = "CpG")
  ) +
    ggplot2::geom_point() +
    ggplot2::ylab("methylation (%)\napparent after quantification") +
    ggplot2::xlab("actual methylation (%)") +
    ggplot2::labs(
      title = "Test-Locus",
      subtitle = "CpG: Test-CpG"
    ) +
    ggplot2::geom_text(
      data = data.frame(),
      ggplot2::aes(x = -Inf,
                   y = c(max(data$CpG),
                         0.95 * max(data$CpG)),
                   hjust = 0, vjust = 1),
      label = lb1,
      parse = F
    )

  # create plot
  outplot <- rBiasCorrection::calibration_plot(
    plotlist = p,
    coef_hyper = coef_hyper,
    coef_cubic = coef_cubic,
    plot_textsize = plot_textsize,
    minmax = FALSE
  )

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = plot_height,
    width = plot_width
  )
}

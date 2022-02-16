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

#' @title calibration_plot helper function
#'
#' @description Internal function to carry out the plotting of the
#'   calibrations curves.
#'
#' @param plotlist A ggplot object containing a regression plot without
#'   regression curves (output of \code{regression_utility()}).
#' @param coef_hyper A list containing the regression parameters of
#'   the hyperbolic regression equation.
#' @param coef_cubic A list containing the regression parameters of
#'   the cubic regression equation.
#' @inheritParams regression_utility
#' @inheritParams createbarerrorplots
#'
#' @return The function returns a list containing calibration plots.
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
          color = "hyperbolic"
        ),
        size = 1.08
      ) +
      ggplot2::stat_function(
        fun = cubic_eq,
        args = list(a = coef_cubic$a,
                    b = coef_cubic$b,
                    c = coef_cubic$c,
                    d = coef_cubic$d),
        geom = "line",
        ggplot2::aes(
          color = "cubic"
        ),
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
          color = "hyperbolic"
        ),
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
        geom = "line",
        ggplot2::aes(
          color = "cubic"
        ),
        size = 1.08
      )
  }

  outplot <- p +
    ggplot2::geom_abline(
      mapping = ggplot2::aes(
        intercept = 0,
        slope = 1,
        color = "unbiased"
      ),
      linetype = "dashed",
      size = 1.08,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      color = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(
      values = c("#E64B35FF",
                 "#4DBBD5FF",
                 "grey47"
                 ),
      labels = c("cubic",
                 "hyperbolic",
                 "unbiased")
    ) +
    #% scale_colour_manual("Regression:",
    #%                     values = c(Cubic = "indianred1",
    #%                                Hyperbolic = "mediumspringgreen",
    #%                                unbiased = "lightblue")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(size = plot_textsize)
    )

  if (isTRUE(getOption("ggpubr.exists"))) {
    outplot <- outplot + ggpubr::theme_pubr()
  }

  return(outplot)
}

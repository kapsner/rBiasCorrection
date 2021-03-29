# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2021 Lorenz Kapsner
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
#' @return The function creates an example plot and stores
#'   on the local filesystem.
#'
#' @examples
#' gdat <- rBiasCorrection::example._plot.df_agg
#'
#' coef_h <- rBiasCorrection::example._plot_coef_h
#' coef_c <- rBiasCorrection::example._plot_coef_c
#'
#' create_exampleplot(
#'   data = gdat,
#'   coef_hyper = coef_h,
#'   coef_cubic = coef_c,
#'   plot_height = 5,
#'   plot_width = 7.5,
#'   plot_textsize = 1,
#'   filename = paste0(tempdir(), "/exampleplot.png")
#' )
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

  data <- data[
    , ("ymin") := get("CpG") - get("sd")
  ][
    , ("ymax") := get("CpG") + get("sd")
  ]

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
      title = "Test Locus",
      subtitle = "CpG: Test CpG"
    ) +
    ggplot2::geom_text(
      data = data.frame(),
      ggplot2::aes(x = -Inf,
                   y = c(100, 0.95 * 100),
                   hjust = 0, vjust = 1),
      label = lb1,
      parse = F
    ) + ggplot2::geom_pointrange(
      ggplot2::aes_string(
        ymin = "ymin",
        ymax = "ymax"
      ),
      fatten = 1
    )

  # create plot
  outplot <- calibration_plot(
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
    width = plot_width,
    dpi = 600
  )
}

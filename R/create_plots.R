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

create_plots <- function(plotlist,
                         f,
                         vec_cal,
                         result_list,
                         filename,
                         logfilename,
                         mode = NULL,
                         minmax,
                         plot_height = 5,
                         plot_width = 7.5,
                         plot_textsize = 16) {

  # hyperbolic parameters
  coef_h <- result_list[[vec_cal[f]]][["Coef_hyper"]]
  # cubic parameters
  coef_c <- result_list[[vec_cal[f]]][["Coef_cubic"]]

  if (isFALSE(minmax)) {

    # create messages
    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", a =", coef_h$a,
                  ", b =", coef_h$b,
                  ", d =", coef_h$d)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Cubic: Using a =", coef_c$a,
                  ", b =", coef_c$b,
                  ", c =", coef_c$c,
                  ", d =", coef_c$d)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    outplot <- calibration_plot(
      plotlist = plotlist,
      coef_hyper = coef_h,
      coef_cubic = coef_c,
      plot_textsize = plot_textsize,
      minmax = minmax
    )

  } else if (isTRUE(minmax)) {
    # create messages
    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Hyperbolic: Using bias_weight =", coef_h$b,
                  ", y0 =", coef_h$y0,
                  ", y1 =", coef_h$y1,
                  ", m0 =", coef_h$m0,
                  ", m1 =", coef_h$m1)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    message <- paste0("# CpG-site: ", vec_cal[f])
    msg2 <- paste("Cubic: Using a =", coef_c$a,
                  ", b =", coef_c$b,
                  ", y0 =", coef_c$y0,
                  ", y1 =", coef_c$y1,
                  ", m0 =", coef_c$m0,
                  ", m1 = ", coef_c$m1)
    write_log(message = paste(message, msg2, sep = "\n"),
              logfilename = logfilename)

    outplot <- calibration_plot(
      plotlist = plotlist,
      coef_hyper = coef_h,
      coef_cubic = coef_c,
      plot_textsize = plot_textsize,
      minmax = minmax
    )
  }

  ggplot2::ggsave(
    filename = filename,
    plot = outplot,
    device = "png",
    height = plot_height,
    width = plot_width,
    dpi = 600
  )
}

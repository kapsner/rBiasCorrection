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


# implementation of cubic equation
cubic_equation <- function(x, a, b, c, d){
  return((a*I(x^3) + b*I(x^2) + c*x + d))
}

cubic_equationMinMax <- function(x, a, b, y0, y1, m0, m1){
  return((a * I((x - m0)^3) + b * I((x - m0)^2) + (x - m0) * (((y1 - y0) / (m1 - m0)) - a * I((m1 - m0)^2) - b * (m1 - m0)) + y0))
}

# find best parameters for cubic regression
cubic_regression <- function(df_agg, vec, logfilename, minmax = minmax) {
  writeLog_("Entered 'cubic_regression'-Function", logfilename)

  dat <- df_agg

  # true y-values
  true_levels <- dat[,get("true_methylation")]

  if (isFALSE(minmax)){
    writeLog_("'cubic_regression': minmax = FALSE", logfilename)

    #pol_reg <- lm(true_methylation ~ poly(CpG, degree = 3, raw = T), data = dat)
    pol_reg <- stats::lm(CpG ~ true_methylation + I(true_methylation^2) + I(true_methylation^3), data = dat)
    cof <- stats::coefficients(pol_reg)

    # correct values
    fitted_values <- cubic_equation(x = true_levels, a = cof[4], b  = cof[3], c = cof[2], d = cof[1])

  } else if (isTRUE(minmax)){
    writeLog_("'cubic_regression': minmax = TRUE --> WARNING: this is experimental", logfilename)

    # extract parameters of equation
    y0 <- dat[get("true_methylation")==dat[,min(get("true_methylation"))], get("CpG")]
    y1 <- dat[get("true_methylation")==dat[,max(get("true_methylation"))], get("CpG")]
    m0 <- dat[,min(get("true_methylation"))]
    m1 <- dat[,max(get("true_methylation"))]

    # starting values
    st <- data.frame(a = c(-1000, 1000),
                     b = c(-1000, 1000))
    set.seed(1234)
    c <- nls2::nls2(CpG ~ cubic_equationMinMax(true_methylation, a, b, y0, y1, m0, m1), data=dat, start = st)

    # get coefficients
    coe <- stats::coef(c)
    a <- coe[["a"]]
    b <- coe[["b"]]

    fitted_values <- cubic_equationMinMax(true_levels, a, b, y0, y1, m0, m1)

  }

  # fitted values
  dat[, ("fitted") := fitted_values]

  # sum of squares between fitted and measuerd values
  dat[,("CpG_fitted_diff") := get("CpG")-get("fitted")]
  dat[,("squared_error") := I((get("CpG_fitted_diff"))^2)]

  # sum of squared errors = residual sum of squares
  SSE <- as.numeric(dat[,sum(get("squared_error"), na.rm = T)])

  # squared dist to mean
  dat[,("squared_dist_mean") := sdm(get("fitted"))]

  # total sum of squares
  TSS <- as.numeric(dat[,sum(get("squared_dist_mean"), na.rm = T)])


  # sum of squared errors
  outlist <- list("SSE_cubic" = SSE)

  if (isFALSE(minmax)){
    outlist[["Coef_cubic"]] = list("ax3" = unname(cof[4]),
                                   "bx2" = unname(cof[3]),
                                   "cx" = unname(cof[2]),
                                   "d" = unname(cof[1]),
                                   "R2" = 1 - (SSE / TSS))
  } else if (isTRUE(minmax)){

    outlist[["Coef_cubic"]] = list("y0" = y0,
                                   "y1" = y1,
                                   "a" = a,
                                   "b" = b,
                                   "m0" = m0,
                                   "m1" = m1,
                                   "R2" = 1 - (SSE / TSS))
  }
  return(outlist)
}

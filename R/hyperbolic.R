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


# implementation of hyperbolic equation
hyperbolic_equationMinMax <- function(x, b, y0, y1, m0, m1){
  # old equation (16.01.2019) with min = 0, max = 100
  #return((((y1 * b) - y0) * x + 100 * y0) / ((b * x) - x + 100))
  # new equation (17.01.2019) with data-dependent min and max
  return((((b * y1) - y0) * (x - m0) + (m1 - m0) * y0) / ((b - 1) * (x - m0) + (m1 - m0)))
}

hyperbolic_equation <- function(x, a, b, d){
  # new equation (24.07.2019) without min and max
  return(((a * x) + b) / (x + d))
}

# solved hyperbolic equation
hyperbolic_equation_solvedMinMax <- function(y, b, y0, y1, m0, m1){
  # old solved equation
  #return(((100 * y0) - (100 * y)) / ((y * b) - (y1 * b) + y0 - y))
  # new solved equation
  return(((m0 * b * (y - y1)) + (m1 * (y0 - y))) / ((b * (y - y1)) - y + y0))
}

hyperbolic_equation_solved <- function(y, a, b, d){
  # new solved equation (24.07.2019) without min and max
  return((b - (y * d)) / (y - a))
}

# find best parameters for hyperbolic regression
hyperbolic_regression <- function(df_agg, vec, logfilename, minmax){
  writeLog_("Entered 'hyperbolic_regression'-Function", logfilename)

  dat <- df_agg

  # true y-values
  true_levels <- dat[,get("true_methylation")]

  # calculate relative error
  dat[,("CpG_true_diff") := abs(get("CpG")-get("true_methylation"))]
  dat[,("relative_error") := ifelse(get("true_methylation") != 0, (get("CpG_true_diff")/get("true_methylation"))*100, NA)]


  if (isFALSE(minmax)){
    writeLog_("'hyperbolic_regression': minmax = FALSE", logfilename)

    # https://cran.r-project.org/web/packages/nls2/nls2.pdf
    # calculate optimal starting values: (longer runtime; same results?)
    # st <- data.frame(a = seq(-5000, 5000, length.out = 2),
    #                  b = seq(-5000, 5000, length.out = 2),
    #                  d = seq(-5000, 5000, length.out = 2))
    # # let's do a grid search to find the best parameters to start with a maximum of 10000 iterations
    # mod <- nls2::nls2(CpG ~ hyperbolic_equation(true_methylation, a, b, d), data=dat, start = st, algorithm = "brute-force", control = nls.control(maxiter = 10000))
    # c <- nls2::nls2(CpG ~ hyperbolic_equation(true_methylation, a, b, d), data=dat, start = mod)


    # starting values
    st <- data.frame(a = c(-1000, 1000),
                     b = c(-1000, 1000),
                     d = c(-1000, 1000))
    set.seed(1234)
    c <- nls2::nls2(CpG ~ hyperbolic_equation(true_methylation, a, b, d), data=dat, start = st)

    # get coefficients
    coe <- stats::coef(c)
    a <- coe[["a"]]
    b <- coe[["b"]]
    d <- coe[["d"]]

    fitted_values <- hyperbolic_equation(true_levels, a, b, d)

  } else if (isTRUE(minmax)){
    writeLog_("'hyperbolic_regression': minmax = TRUE --> WARNING: this is experimental", logfilename)

    # extract parameters of equation
    y0 <- dat[get("true_methylation")==dat[,min(get("true_methylation"))], get("CpG")]
    y1 <- dat[get("true_methylation")==dat[,max(get("true_methylation"))], get("CpG")]
    m0 <- dat[,min(get("true_methylation"))]
    m1 <- dat[,max(get("true_methylation"))]

    # implementation of optimization function
    fn <- function(bias){
      fitted_vals <- hyperbolic_equationMinMax(true_levels, b = bias, y0 = y0, y1 = y1, m0 = m0, m1 = m1)
      # optimize biasfactor with minimizing sum of squares error
      return(sum(I(dat[,get("CpG")] - fitted_vals)^2))
    }

    # optimization function of built in R -> based on Nelder-Mead
    # by default, optim performs minimization
    # bias_factor <- optim(1, fn, method = "Nelder-Mead")$par
    b <- stats::optim(1, fn, method = "Brent", lower = 0, upper = 50)$par # due to error with Nelder-Mead

    # correct values, based on optimized b
    fitted_values <- hyperbolic_equationMinMax(true_levels, b, y0, y1, m0, m1)
  }

  # the next part is equal for minmax = FALSE and minmax = TRUE

  # fitted values, extrapolated by true methylation and y0 and y1
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
  outlist <- list("Var" = vec,
                  "relative_error" = dat[,mean(get("relative_error"), na.rm = T)],
                  "SSE_hyper" = SSE)

  if (isFALSE(minmax)){
    outlist[["Coef_hyper"]] = list("a" = a,
                                   "b" = b,
                                   "d" = d,
                                   "R2" = 1 - (SSE / TSS))
  } else if (isTRUE(minmax)){

    outlist[["Coef_hyper"]] = list("y0" = y0,
                                   "y1" = y1,
                                   "b" = b,
                                   "m0" = m0,
                                   "m1" = m1,
                                   "R2" = 1 - (SSE / TSS))
  }
  return(outlist)
}

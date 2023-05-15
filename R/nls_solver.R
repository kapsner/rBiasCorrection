# rBiasCorrection: Correct Bias in Quantitative DNA Methylation Analyses.
# Copyright (C) 2019-2023 Lorenz Kapsner
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

nls_solver <- function(
    true_levels,
    target_levels,
    type = c("hyperbolic_eq", "hyperbolic_eq_minmax", "cubic_eq_minmax"),
    logfilename,
    seed,
    ...) {

  type <- match.arg(type)

  if (type %in% c("hyperbolic_eq_minmax", "cubic_eq_minmax")) {
    required_dot_args <- c("y0", "y1", "m0", "m1")
    kwargs <- list(...)
    kwargs_names <- names(kwargs)
    for (da in required_dot_args) {
      if (da %in% kwargs_names) {
        assign(
          x = da,
          value = kwargs[[da]]
        )
      } else {
        stop(paste0("'", da, "' is missing in dot-args ..."))
      }
    }
  }

  formula_common <- "target_levels ~ FUN(x = true_levels, b = b, "
  minmax_common <- "y0 = y0, y1 = y1, m0 = m0, m1 = m1)"

  FUN <- switch( # nolint
    EXPR = type,
    "hyperbolic_eq" = hyperbolic_eq,
    "hyperbolic_eq_minmax" = hyperbolic_eq_minmax,
    "cubic_eq_minmax" = cubic_eq_minmax
  )

  FUN_formula <- switch( # nolint
    EXPR = type,
    "hyperbolic_eq" = paste0(formula_common, "a = a, d = d)"),
    "hyperbolic_eq_minmax" = paste0(formula_common, minmax_common),
    "cubic_eq_minmax" = paste0(formula_common, "a = a, ", minmax_common)
  )
  FUN_formula <- stats::as.formula( # nolint
    object = FUN_formula
  )

  nls_switch <- getOption("rBiasCorrection.nls_implementation")
  nls_options <- c("nls2_paper", "nls2_fast", "minpack.lm")
  if (is.null(nls_switch) || !(nls_switch %in% nls_options)) {
    nls_switch <- "nls2_paper"
  }

  if (nls_switch %in% grep("nls", nls_options, value = TRUE)) {

    # starting values
    if (nls_switch == "nls2_paper") {
      start_vals <- switch( # nolint
        EXPR = type,
        "hyperbolic_eq" = data.frame(a = c(-1000, 1000),
                                     b = c(-1000, 1000),
                                     d = c(-1000, 1000)),
        "hyperbolic_eq_minmax" = data.frame(b = c(-1000, 1000)),
        "cubic_eq_minmax" = data.frame(a = c(-1000, 1000),
                                       b = c(-1000, 1000))
      )
    } else if (nls_switch == "nls2_fast") {
      coef_df <- fast_guess(
        type = type,
        target_levels = target_levels,
        true_levels = true_levels
      )

      start_vals <- switch( # nolint
        EXPR = type,
        "hyperbolic_eq" = data.frame(
          a = c(coef_df[["a"]] * -1, coef_df[["a"]]),
          b = c(coef_df[["b"]] * -1, coef_df[["b"]]),
          d = c(coef_df[["d"]] * -1, coef_df[["d"]])
        ),
        "hyperbolic_eq_minmax" = data.frame(
          b = c(coef_df[["b"]] * -1, coef_df[["b"]])
        ),
        "cubic_eq_minmax" = data.frame(
          a = c(coef_df[["a"]] * -1, coef_df[["a"]]),
          b = c(coef_df[["b"]] * -1, coef_df[["b"]])
        )
      )
    }

    c <- tryCatch({
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
      set.seed(seed)
      ret <- nls2::nls2(
        formula = FUN_formula,
        start = start_vals,
        control = stats::nls.control(maxiter = 50)
      )
      ret

    }, error = function(e) {
      # if convergence fails
      write_log(message = e,
                logfilename = logfilename)
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
      set.seed(seed)
      mod <- nls2::nls2(
        formula = FUN_formula,
        start = start_vals,
        algorithm = "brute-force",
        control = stats::nls.control(maxiter = 1e5)
      )

      suppressWarnings(RNGkind(sample.kind = "Rounding"))
      set.seed(seed)
      ret <- nls2::nls2(
        formula = FUN_formula,
        start = mod,
        algorithm = "brute-force",
        control = stats::nls.control(maxiter = 1e3)
      )
      ret
    }, finally = function(f) {
      return(ret)
    })

  } else if (nls_switch == "minpack.lm") {

    start_vals <- fast_guess(
      type = type,
      target_levels = target_levels,
      true_levels = true_levels
    )

    c <- minpack.lm::nlsLM(
      formula = FUN_formula,
      start = start_vals,
      algorithm = "LM",
      control = stats::nls.control(maxiter = 50)
    )
  } else {
    stop("Not implemented.")
  }

  return(c)
}

fast_guess <- function(type, target_levels, true_levels) {
  if (grepl("hyperbolic", type)) {
    d <- 0.001
    guess_mod <- guess_nls_start_linear(
      target_levels = target_levels,
      true_levels = true_levels
    )
  } else if (type == "cubic_eq_minmax") {
    guess_mod <- cubic_fitter(
      target_levels = target_levels,
      true_levels = true_levels
    )
  }

  mod_coef <- stats::coef(guess_mod)

  if (type == "hyperbolic_eq") {
    # a, b, d
    # b is intercept in hyperbolic_eq
    coef_ret <- c(mod_coef[2], mod_coef[1], d)
    names(coef_ret) <- c("a", "b", "d")
  } else if (type == "hyperbolic_eq_minmax") {
    # b is beta
    coef_ret <- mod_coef[2]
    names(coef_ret) <- "b"
  } else if (type == "cubic_eq_minmax") {
    # a = x3, b = x2
    coef_ret <- c(mod_coef[4], mod_coef[3])
    names(coef_ret) <- c("a", "b")
  }
  return(coef_ret)
}

guess_nls_start_linear <- function(target_levels, true_levels) {
  lin_mod <- stats::lm(target_levels ~ true_levels)
  return(lin_mod)
}

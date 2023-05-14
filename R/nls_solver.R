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
    "hyperbolic_eq" = as.formula(
      paste0(formula_common, "a = a, d = d)")
    ),
    "hyperbolic_eq_minmax" = as.formula(
      paste0(formula_common, minmax_common)
    ),
    "cubic_eq_minmax" = as.formula(
      paste0(formula_common, "a = a, ", minmax_common)
    )
  )

  nls_switch <- getOption("rBiasCorrection.nls_implementation")
  if (is.null(nls_switch) || nls_switch == "") {
    nls_switch <- "nls2"
  }

  if (nls_switch == "nls2") {

    # starting values
    start <- switch( # nolint
      EXPR = type,
      "hyperbolic_eq" = data.frame(a = c(-1000, 1000),
                                   b = c(-1000, 1000),
                                   d = c(-1000, 1000)),
      "hyperbolic_eq_minmax" = data.frame(b = c(-1000, 1000)),
      "cubic_eq_minmax" = data.frame(a = c(-1000, 1000),
                                     b = c(-1000, 1000))
    )

    c <- tryCatch({
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
      set.seed(seed)
      ret <- nls2::nls2(
        formula = FUN_formula,
        start = start,
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
        start = start,
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

    lin_mod <- lm(target_levels ~ true_levels)
    lin_mod_coef <- stats::coef(lin_mod)
    coef_df <- c(lin_mod_coef[2], lin_mod_coef[1], 0.001)
    names(coef_df) <- c("a", "b", "d")

    start <- switch( # nolint
      EXPR = type,
      "hyperbolic_eq" = coef_df,
      "hyperbolic_eq_minmax" = coef_df["b"],
      "cubic_eq_minmax" = coef_df[c("a", "b")]
    )

    c <- minpack.lm::nlsLM(
      formula = FUN_formula,
      start = start,
      algorithm = "LM",
      control = stats::nls.control(maxiter = 50)
    )
  } else {
    stop("Not implemented.")
  }

  return(c)
}

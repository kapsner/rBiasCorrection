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

  # starting values
  # TODO implement logic to better guess of initial starting values
  start <- switch( # nolint
    EXPR = type,
    "hyperbolic_eq" = data.frame(a = c(-1000, 1000),
                                 b = c(-1000, 1000),
                                 d = c(-1000, 1000)),
    "hyperbolic_eq_minmax" = data.frame(b = c(-1000, 1000)),
    "cubic_eq_minmax" = data.frame(a = c(-1000, 1000),
                                   b = c(-1000, 1000))
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
    st <- sapply(
      X = colnames(start),
      FUN = function(x) {
        start[1, x]
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    c <- minpack.lm::nlsLM(
      formula = FUN_formula,
      start = st,
      algorithm = "LM",
      control = stats::nls.control(maxiter = 50)
    )
  } else {
    stop("Not implemented.")
  }

  return(c)
}

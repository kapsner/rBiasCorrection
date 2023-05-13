nls_solver <- function(
    true_levels,
    target_levels,
    start,
    type = c("hyperbolic_eq", "hyperbolic_eq_minmax"),
    seed,
    ...) {

  type <- match.arg(type)

  formula_common <- "target_levels ~ FUN(x = true_levels, b = b, "

  FUN <- switch(
    EXPR = type,
    "hyperbolic_eq" = hyperbolic_eq,
    "hyperbolic_eq_minmax" = hyperbolic_eq_minmax
  )

  FUN_formula <- switch(
    EXPR = type,
    "hyperbolic_eq" = as.formula(
      paste0(formula_common, "a = a, d = d)")
    ),
    "hyperbolic_eq_minmax" = as.formula(
      paste0(formula_common, "y0 = y0, y1 = y1, m0 = m0, m1 = m1)")
    )
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
      start = st,
      algorithm = "brute-force",
      control = stats::nls.control(maxiter = 1e5)
    )

    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    set.seed(seed)
    ret <- nls2::nls2(
      formula = FUN_formula,
      data = dat,
      start = mod,
      algorithm = "brute-force",
      control = stats::nls.control(maxiter = 1e3)
    )
    ret
  }, finally = function(f) {
    return(ret)
  })
  return(c)
}

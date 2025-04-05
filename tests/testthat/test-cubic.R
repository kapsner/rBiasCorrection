prefix <- tempdir()
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")


library(data.table)
# calibration data
cal_type_1 <- data.table::fread("./testdata/cal_type_1.csv")
# cal_type_1 <- data.table::fread(paste0(prefix, "/testdata/cal_type_1.csv"))
cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])

test_that(
  desc = "test functioning of cubic regression function - nls2",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    h1 <- cubic_regression(df_agg = df_agg,
                           vec = colnames(cal_type_1)[2],
                           logfilename,
                           minmax = TRUE,
                           seed = 1234)
    h1_rounded <- list(
      SSE_cubic = h1$SSE_cubic,
      Coef_cubic = lapply(h1$Coef_cubic, round, digits = 2)
    )
    expect_snapshot(
      x = h1_rounded,
      cran = FALSE,
      error = FALSE
    )

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })


test_that(
  desc = "test functioning of cubic regression function - minpack.lm",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    options(rBiasCorrection.nls_implementation = "LM")
    h2 <- cubic_regression(df_agg = df_agg,
                           vec = colnames(cal_type_1)[2],
                           logfilename,
                           minmax = TRUE,
                           seed = 1234)

    h2_rounded <- list(
      SSE_cubic = h2$SSE_cubic,
      Coef_cubic = lapply(h2$Coef_cubic, round, digits = 2)
    )

    expect_snapshot(
      x = h2_rounded,
      cran = FALSE,
      error = FALSE
    )
    # reset
    options(rBiasCorrection.nls_implementation = "GN.paper")

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

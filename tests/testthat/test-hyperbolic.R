prefix <- tempdir()
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")


library(data.table)
# calibration data
cal_type_1 <- data.table::fread("./testdata/cal_type_1.csv")
cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])

test_that(
  desc = "test functioning of hyperbolic regression function - nls2",
  code = {

    local_edition(3)
    #"skip_on_cran()

    h1 <- hyperbolic_regression(df_agg = df_agg,
                                vec = colnames(cal_type_1)[2],
                                logfilename,
                                minmax = TRUE,
                                seed = 1234)
    expect_snapshot_value(
      x = h1,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })


test_that(
  desc = "test functioning of hyperbolic regression function - minpack.lm",
  code = {

    local_edition(3)
    #"skip_on_cran()

    options(rBiasCorrection.nls_implementation = "LM")
    h2 <- hyperbolic_regression(df_agg = df_agg,
                                vec = colnames(cal_type_1)[2],
                                logfilename,
                                minmax = TRUE,
                                seed = 1234)
    expect_snapshot_value(
      x = h2,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    options(rBiasCorrection.nls_implementation = "GN.paper")

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

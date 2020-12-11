context("test hyperbolic")

prefix <- tempdir()
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")


library(data.table)

test_that(
  desc = "test functioning of hyperbolic regression function",
  code = {

    local_edition(3)
    #"skip_on_cran()

    # calibration data
    cal_type_1 <- fread("./testdata/cal_type_1.csv")
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
    df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])

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
    # nolint start
    # expect_error({
    #   expect_known_hash(h1, "3578c7d484")
    #   expect_known_hash(h1, "bce2d004e7")
    # }, class = "error", regexp = "3578c7d484|bce2d004e7") # bce2d004e7
    # nolint end

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

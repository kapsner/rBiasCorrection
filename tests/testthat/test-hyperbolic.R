context("test hyperbolic")

prefix <- "./"
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")


library(data.table)

test_that(
  desc = "test functioning of hyperbolic regression function",
  code = {

    skip_on_cran()

    # calibration data
    cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
    df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])

    h1 <- hyperbolic_regression(df_agg = df_agg,
                                vec = colnames(cal_type_1)[2],
                                logfilename,
                                minmax = TRUE,
                                seed = 1234)
    expect_known_hash(h1, "bce2d004e7") # 3578c7d484
  })

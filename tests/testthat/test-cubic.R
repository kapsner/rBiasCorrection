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
  desc = "test functioning of cubic regression function - nls2",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    h1 <- cubic_regression(df_agg = df_agg,
                           vec = colnames(cal_type_1)[2],
                           logfilename,
                           minmax = TRUE,
                           seed = 1234)
    expect_snapshot(
      x = h1,
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
    expect_snapshot(
      x = h2,
      cran = FALSE,
      error = FALSE
    )
    # reset
    options(rBiasCorrection.nls_implementation = "GN.paper")

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

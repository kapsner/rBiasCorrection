prefix <- tempdir()
# prefix <- "tests/testthat/" # nolint

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")


library(data.table)

test_that(
  desc = "test functioning of aggregated function",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    # calibration data
    cal_type_1 <- rBiasCorrection::example.data_calibration$dat
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
    df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])
    expect_snapshot(
      x = round(table_prep(df_agg), 1),
      cran = FALSE,
      error = FALSE
    )

    # experimental data
    exp_type_1 <- rBiasCorrection::example.data_experimental$dat
    exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
    df_agg <- create_agg_df_exp(exp_type_1, colnames(exp_type_1)[2], type = 1)
    expect_snapshot(
      x = round(table_prep(df_agg), 1),
      cran = FALSE,
      error = FALSE
    )

    exp_type_2 <- fread("./testdata/exp_type_2.csv")
    exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
    df_agg <- create_agg_df_exp(exp_type_2, colnames(exp_type_2)[2], type = 2)
    expect_snapshot(
      x = round(table_prep(df_agg), 2),
      cran = FALSE,
      error = FALSE
    )

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

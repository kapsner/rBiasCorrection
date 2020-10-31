context("create df_agg")

prefix <- tempdir()
# prefix <- "tests/testthat/" # nolint

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")


library(data.table)

test_that(
  desc = "test functioning of aggregated function",
  code = {
  # calibration data
  cal_type_1 <- rBiasCorrection::example.data_calibration$dat
  cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
  df_agg <- create_agg_df(cal_type_1, colnames(cal_type_1)[2])
  expect_error({
    expect_known_hash(df_agg, "34d2e049fb")
    expect_known_hash(df_agg, "c16660dcd8")
  }, class = "error", regexp = "34d2e049fb|c16660dcd8")  # 6aa2d6fc51

  # experimental data
  exp_type_1 <- rBiasCorrection::example.data_experimental$dat
  exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
  df_agg <- create_agg_df_exp(exp_type_1, colnames(exp_type_1)[2], type = 1)
  expect_error({
    expect_known_hash(df_agg, "5a885b7e0a")
    expect_known_hash(df_agg, "4a4ed97b38")
  }, class = "error", regexp = "5a885b7e0a|4a4ed97b38")

  exp_type_2 <- fread("./testdata/exp_type_2.csv")
  exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
  df_agg <- create_agg_df_exp(exp_type_2, colnames(exp_type_2)[2], type = 2)
  expect_error({
    expect_known_hash(df_agg, "d042a19f13")
    expect_known_hash(df_agg, "46d527e146")
  }, class = "error", regexp = "d042a19f13|46d527e146")

  expect_true(file.remove(paste0(prefix, "/log.txt")))
})

prefix <- tempdir()
# prefix <- "tests/testthat/" # nolint

library(data.table)

test_that(
  desc = "correct functioning of BiasCorrection, data type 1",
  code = {

  plotdir <- paste0(prefix, "/plotdir/")
  csvdir <- paste0(prefix, "/csvdir/")

  b1 <- biascorrection(
    experimental = "./testdata/exp_type_1.csv",
    calibration = "./testdata/cal_type_1.csv",
    samplelocusname = "Testlocus",
    minmax = FALSE,
    correct_method = "best",
    selection_method = "SSE",
    type = 1,
    plotdir = plotdir,
    csvdir = csvdir,
    logfilename = paste0(prefix, "/log.txt"),
    seed = 1234,
    parallel = ifelse(
      tolower(Sys.info()["sysname"]) == "darwin",
      FALSE,
      TRUE
    )
  )
  expect_type(
    b1,
    "list"
  )
  expect_length(list.files(plotdir), 50)
  expect_length(list.files(csvdir), 13)

  options(rBiasCorrection.nls_implementation = "nls2_fast")
  b2 <- biascorrection(
    experimental = "./testdata/exp_type_1.csv",
    calibration = "./testdata/cal_type_1.csv",
    samplelocusname = "Testlocus",
    minmax = FALSE,
    correct_method = "best",
    selection_method = "SSE",
    type = 1,
    plotdir = plotdir,
    csvdir = csvdir,
    logfilename = paste0(prefix, "/log.txt"),
    seed = 1234,
    parallel = ifelse(
      tolower(Sys.info()["sysname"]) == "darwin",
      FALSE,
      TRUE
    )
  )

  options(rBiasCorrection.nls_implementation = "minpack.lm")
  b3 <- biascorrection(
    experimental = "./testdata/exp_type_1.csv",
    calibration = "./testdata/cal_type_1.csv",
    samplelocusname = "Testlocus",
    minmax = FALSE,
    correct_method = "best",
    selection_method = "SSE",
    type = 1,
    plotdir = plotdir,
    csvdir = csvdir,
    logfilename = paste0(prefix, "/log.txt"),
    seed = 1234,
    parallel = ifelse(
      tolower(Sys.info()["sysname"]) == "darwin",
      FALSE,
      TRUE
    )
  )

  expect_equal(
    object = b3$final_results,
    expected = b2$final_results
  )
  options(rBiasCorrection.nls_implementation = "nls2_paper")

  # cleanup
  expect_silent(clean_up(plotdir = plotdir,
                         csvdir = csvdir))
  expect_true(file.remove(paste0(prefix, "/log.txt")))
})

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
    correct_method = "hyperbolic",
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
  expect_length(list.files(plotdir), 50)
  expect_length(list.files(csvdir), 13)

  options(rBiasCorrection.nls_implementation = "GN.guess")
  b2 <- biascorrection(
    experimental = "./testdata/exp_type_1.csv",
    calibration = "./testdata/cal_type_1.csv",
    samplelocusname = "Testlocus",
    minmax = FALSE,
    correct_method = "hyperbolic",
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
    object = b1$final_results,
    expected = b2$final_results,
    tolerance = 1e-3
  )

  options(rBiasCorrection.nls_implementation = "LM")
  b3 <- biascorrection(
    experimental = "./testdata/exp_type_1.csv",
    calibration = "./testdata/cal_type_1.csv",
    samplelocusname = "Testlocus",
    minmax = FALSE,
    correct_method = "hyperbolic",
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
    object = b1$final_results,
    expected = b3$final_results,
    tolerance = 2e-1
  )
  options(rBiasCorrection.nls_implementation = "GN.paper")

  # cleanup
  expect_silent(clean_up(plotdir = plotdir,
                         csvdir = csvdir))
  expect_true(file.remove(paste0(prefix, "/log.txt")))
})

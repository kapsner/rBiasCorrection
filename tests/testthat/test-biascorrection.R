context("test BiasCorrection function")

prefix <- "./"
#" prefix <- "tests/testthat/"

library(data.table)

test_that(
  desc = "correct functioning of BiasCorrection, data type 1",
  code = {

  plotdir <- paste0(prefix, "plotdir")
  csvdir <- paste0(prefix, "csvdir")

  expect_true(
    biascorrection(
      experimental = paste0(prefix, "testdata/exp_type_1.csv"),
      calibration = paste0(prefix, "testdata/cal_type_1.csv"),
      samplelocusname = "Testlocus",
      minmax = FALSE,
      correct_method = "best",
      selection_method = "SSE",
      type = 1,
      plotdir = plotdir,
      csvdir = csvdir,
      logfilename = paste0(prefix, "log.txt"),
      seed = 1234
    )
  )
  expect_length(list.files(plotdir), 50)
  expect_length(list.files(csvdir), 11)


  # cleanup
  expect_silent(clean_up(plotdir = plotdir,
                         csvdir = csvdir))
  expect_true(file.remove(paste0(prefix, "log.txt")))
})

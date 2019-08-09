context("test BiasCorrection function")

prefix <- "./"
#prefix <- "tests/testthat/"

library(data.table)

test_that("correct functioning of BiasCorrection, data type 1",{

  plotdir <- paste0(prefix, "plotdir")
  csvdir <- paste0(prefix, "csvdir")

  expect_true(BiasCorrection(experimental = paste0(prefix, "testdata/exp_type_1.csv"),
                             calibration = paste0(prefix, "testdata/cal_type_1.csv"),
                             samplelocusname = "Testlocus",
                             minmax = FALSE,
                             correct_method = "best",
                             selection_method = "SSE",
                             type = 1,
                             plotdir = plotdir,
                             csvdir = csvdir,
                             logfilename = paste0(prefix, "log.txt")))
  expect_length(list.files(plotdir), 50)
  expect_length(list.files(csvdir), 11)


  # cleanup
  expect_silent(cleanUp_(plotdir = plotdir,
                         csvdir = csvdir))
  expect_true(file.remove(paste0(prefix, "log.txt")))
})




# debug all
# prefix <- "tests/testthat/"
# plotdir <- paste0(prefix, "plotdir")
# csvdir <- paste0(prefix, "csvdir")
# experimental = paste0(prefix, "testdata/exp_type_1.csv")
# calibration = paste0(prefix, "testdata/cal_type_1.csv")
# samplelocusname = "Testlocus"
# minmax = TRUE
# type = 1
# logfilename = paste0(prefix, "log.txt")

# debug solvingEquations_
# datatable = rv$fileimportCal
# regmethod = rv$choices_list
# type = 1
# rv = rv
# mode = "corrected"
# logfilename = logfilename
# minmax = rv$minmax

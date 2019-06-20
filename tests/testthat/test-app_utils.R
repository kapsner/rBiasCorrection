context("test app utils")

prefix <- "./"

library(data.table)


test_that("correct functioning app utils",{

  expect_true(onStart(plotdir = paste0(prefix, "plotdir"),
                        csvdir = paste0(prefix, "csvdir"),
                        logfilename = paste0(prefix, "log.txt")))


  # cleanup
  expect_silent(cleanUp(plotdir = paste0(prefix, "plotdir"),
                        csvdir = paste0(prefix, "csvdir")))
  expect_true(file.remove(paste0(prefix, "log.txt")))
})

context("test app utils")

prefix <- "./"

library(data.table)


test_that(
  desc = "correct functioning app utils",
  code = {

    expect_silent(
      on_start(plotdir = paste0(prefix, "plotdir"),
               csvdir = paste0(prefix, "csvdir"),
               logfilename = paste0(prefix, "log.txt"),
               parallel = TRUE)
    )


    # cleanup
    expect_silent(clean_up(plotdir = paste0(prefix, "plotdir"),
                           csvdir = paste0(prefix, "csvdir")))
    expect_true(file.remove(paste0(prefix, "log.txt")))
  })

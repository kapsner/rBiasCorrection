prefix <- tempdir()

library(data.table)


test_that(
  desc = "correct functioning app utils",
  code = {

    plotdir <- paste0(prefix, "/plotdir/")
    csvdir <- paste0(prefix, "/csvdir/")

    vitual_onstart <- on_start(plotdir = plotdir,
                               csvdir = csvdir,
                               logfilename = paste0(prefix, "/log.txt"),
                               parallel = FALSE)
    expect_type(vitual_onstart, type = "list")


    # cleanup
    expect_silent(clean_up(plotdir = plotdir,
                           csvdir = csvdir))
    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

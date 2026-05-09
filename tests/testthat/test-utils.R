prefix <- tempdir()

library(data.table)


test_that(
  desc = "correct functioning app utils",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    plotdir <- file.path(prefix, "plotdir/")
    csvdir <- file.path(prefix, "csvdir/")

    vitual_onstart <- on_start(plotdir = plotdir,
                               csvdir = csvdir,
                               logfilename = file.path(prefix, "log.txt"),
                               parallel = FALSE)
    expect_type(vitual_onstart, type = "list")


    # cleanup
    expect_silent(clean_up(plotdir = plotdir,
                           csvdir = csvdir))
    expect_true(file.remove(file.path(prefix, "log.txt")))
  })

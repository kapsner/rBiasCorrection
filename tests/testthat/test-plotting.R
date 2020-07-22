context("test plotting")

prefix <- "./"
# prefix <- "tests/testthat/" # nolint

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that(
  desc = "plotting_utility",
  code = {

    #"skip_on_cran()

    rv$minmax <- FALSE
    rv$sample_locus_name <- "Test"
    rv$seed <- 1234

    # create plotdir
    plotdir <- paste0(prefix, "plotdir/")
    csvdir <- paste0(prefix, "csvdir/")

    on_start(plotdir = plotdir,
             csvdir = csvdir,
             logfilename = paste0(prefix, "log.txt"),
             parallel = FALSE)

    # experimental data
    exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
    rv$fileimport_experimental <- clean_dt(exp_type_1,
                                           "experimental",
                                           1,
                                           logfilename)[["dat"]]

    # calibration data
    cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)
    rv$fileimport_calibration <- cal_type_1[["dat"]]
    rv$vec_cal <- cal_type_1[["vec_cal"]]

    # reconstruct parts from app_plottingUtility.R
    regression_results <- regression_utility(
      data = rv$fileimport_calibration,
      samplelocusname = "Testlocus",
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

    # create calibration plots
    virtual_list <- plotting_utility(
      data = rv$fileimport_calibration,
      plotlist_reg = plotlist_reg,
      type = 1,
      samplelocusname = rv$sample_locus_name,
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      plotdir = plotdir,
      logfilename = logfilename,
      minmax = rv$minmax,
      plot_height = 5,
      plot_width = 7.5,
      plot_textsize = 1
    )

    expect_type(virtual_list, "list")


    expect_length(list.files(plotdir), 10)

    # cleanup
    expect_silent(clean_up(plotdir = plotdir,
                           csvdir = csvdir))
})

test_that(
  desc = "createbarerrorplots",
  code = {

    #"skip_on_cran()

    rv$minmax <- FALSE
    rv$sample_locus_name <- "Test"
    rv$seed <- 1234

    # create plotdir
    plotdir <- paste0(prefix, "plotdir/")
    csvdir <- paste0(prefix, "csvdir/")

    on_start(plotdir = plotdir,
             csvdir = csvdir,
             logfilename = paste0(prefix, "log.txt"),
             parallel = FALSE)

    # experimental data
    exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
    rv$fileimport_experimental <- clean_dt(exp_type_1,
                                           "experimental",
                                           1,
                                           logfilename)[["dat"]]

    # calibration data
    cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)
    rv$fileimport_calibration <- cal_type_1[["dat"]]
    rv$vec_cal <- cal_type_1[["vec_cal"]]

    # reconstruct parts from app_plottingUtility.R
    regression_results <- regression_utility(
      data = rv$fileimport_calibration,
      samplelocusname = "Testlocus",
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )
    rv$result_list <- regression_results[["result_list"]]

    # get regression statistics
    rv$reg_stats <- statistics_list(
      rv$result_list,
      minmax = rv$minmax
    )

     # select the better model based on the sum of squared errrors ("SSE")
     rv$choices_list <- better_model(
       statstable_pre = rv$reg_stats,
       selection_method = "SSE"
     )

     # correct calibration data (to show corrected calibration curves)
     solved_eq_h <- solving_equations(datatable = rv$fileimport_calibration,
                                      regmethod = rv$choices_list,
                                      type = 1,
                                      rv = rv,
                                      mode = "corrected",
                                      logfilename = logfilename,
                                      minmax = rv$minmax)
     rv$fileimport_cal_corrected_h <- solved_eq_h$results
     colnames(rv$fileimport_cal_corrected_h) <- colnames(
       rv$fileimport_calibration
     )

     # calculate new calibration curves from corrected calibration data
     regression_results <- regression_utility(
       data = rv$fileimport_cal_corrected_h,
       samplelocusname = rv$sample_locus_name,
       rv = rv,
       mode = "corrected",
       logfilename = logfilename,
       minmax = rv$minmax,
       seed = rv$seed
     )
     rv$result_list_hyperbolic <- regression_results$result_list


     # save regression statistics to reactive value
     rv$reg_stats_corrected_h <- statistics_list(
       resultlist = rv$result_list_hyperbolic,
       minmax = rv$minmax
     )

     virtual_list <- createbarerrorplots(
       statstable_pre = rv$reg_stats,
       statstable_post = rv$reg_stats_corrected_h,
       rv = rv,
       type = 1,
       locus_id = NULL,
       plotdir = plotdir,
       logfilename = logfilename,
       mode = "corrected_h",
       plot_height = 5,
       plot_width = 7.5,
       plot_textsize = 1
     )

    expect_type(virtual_list, "list")


    expect_length(list.files(plotdir), 10)

    # cleanup
    expect_silent(clean_up(plotdir = plotdir,
                           csvdir = csvdir))
    expect_true(file.remove(paste0(prefix, "log.txt")))
  })

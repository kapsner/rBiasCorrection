prefix <- tempdir()
# prefix <- "tests/testthat/" # nolint

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that(
  desc = "algorithm test, type 1, minmax = TRUE",
  code = {

    local_edition(3)
    local_reproducible_output(rstudio = TRUE)

    suppressWarnings(future::plan("multisession"))

    #"skip_on_cran()

    rv$minmax <- TRUE
    rv$sample_locus_name <- "Test"
    rv$seed <- 1234

    # experimental data
    exp_type_1 <- rBiasCorrection::example.data_experimental$dat
    rv$fileimport_experimental <- clean_dt(exp_type_1,
                                           "experimental",
                                           1,
                                           logfilename)[["dat"]]

    # calibration data
    cal_type_1 <- rBiasCorrection::example.data_calibration$dat
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)
    rv$fileimport_calibration <- cal_type_1[["dat"]]
    rv$vec_cal <- cal_type_1[["vec_cal"]]

    # some tests
    expect_length(rv$vec_cal, 10)
    expect_type(rv$vec_cal, "character")


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

    options(rBiasCorrection.nls_implementation = "GN.guess")
    regression_results_fast <- regression_utility(
      data = rv$fileimport_calibration,
      samplelocusname = "Testlocus",
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )

    expect_equal(
      object = regression_results$result_list,
      expected =  regression_results_fast$result_list
    )

    options(rBiasCorrection.nls_implementation = "LM")
    regression_results_minpack <- regression_utility(
      data = rv$fileimport_calibration,
      samplelocusname = "Testlocus",
      locus_id = NULL,
      rv = rv,
      mode = NULL,
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )

    expect_equal(
      object = regression_results$result_list,
      expected =  regression_results_minpack$result_list
    )

    options(rBiasCorrection.nls_implementation = "GN.paper")


    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list <- testhelper_apply_robust_results_list(
      regression_results[["result_list"]]
    )

    regression_results2 <- regression_type1(rv$fileimport_calibration,
                                            rv$vec_cal,
                                            mode = NULL,
                                            logfilename,
                                            minmax = rv$minmax,
                                            locus_id = NULL,
                                            locusname = rv$sample_locus_name,
                                            seed = rv$seed)

    # save regression statistics to reactive value
    rv$reg_stats <- statistics_list(rv$result_list, minmax = rv$minmax)

    # some tests
    expect_snapshot(
      x = round(table_prep(rv$reg_stats), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_equal(
        regression_results$result_list,
        regression_results2$result_list
    )

    # calculate final results
    rv$choices_list <- better_model(statstable_pre = rv$reg_stats,
                                    selection_method = "SSE")
    solved_eq <- solving_equations(rv$fileimport_experimental,
                                   rv$choices_list,
                                   type = 1,
                                   rv = rv,
                                   logfilename = logfilename,
                                   minmax = rv$minmax)
    rv$final_results <- solved_eq[["results"]]
    rv$substitutions <- solved_eq[["substitutions"]]

    # Calibration Data (to show corrected calibration curves)
    solved_eq2 <- solving_equations(rv$fileimport_calibration,
                                    rv$choices_list,
                                    type = 1,
                                    rv = rv,
                                    mode = "corrected",
                                    logfilename = logfilename,
                                    minmax = rv$minmax)
    rv$fileimport_cal_corrected <- solved_eq2[["results"]]
    colnames(rv$fileimport_cal_corrected) <-
      colnames(rv$fileimport_calibration)

    # some tests
    expect_snapshot(
      x = round(table_prep(rv$final_results), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(rv$substitutions), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(solved_eq2[["results"]]), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(solved_eq2[["substitutions"]]), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(rv$fileimport_cal_corrected), 2),
      cran = FALSE,
      error = FALSE
    )



    # hyperbolic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = FALSE
    ][
      , ("better_model") := 0
    ]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_h <- solving_equations(rv$fileimport_calibration,
                                     rv$choices_list,
                                     type = 1,
                                     rv = rv,
                                     mode = "corrected",
                                     logfilename = logfilename,
                                     minmax = rv$minmax)
    rv$fileimport_cal_corrected_h <- solved_eq_h[["results"]]
    colnames(rv$fileimport_cal_corrected_h) <- colnames(
      rv$fileimport_calibration
    )
    rv$substitutions_corrected_h <- solved_eq_h[["substitutions"]]

    expect_snapshot(
      x = round(table_prep(rv$fileimport_cal_corrected_h), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(rv$substitutions_corrected_h), 2),
      cran = FALSE,
      error = FALSE
    )

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      rv$fileimport_cal_corrected_h,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list_hyperbolic <- testhelper_apply_robust_results_list(
      regression_results[["result_list"]]
    )
    # save regression statistics to reactive value
    rv$reg_stats_corrected_h <- statistics_list(rv$result_list_hyperbolic,
                                                minmax = rv$minmax)

    expect_snapshot(
      x = round(table_prep(rv$reg_stats_corrected_h), 2),
      cran = FALSE,
      error = FALSE
    )


    # cubic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = FALSE
    ][
      , ("better_model") := 1
    ]

    # correct calibration data (to show corrected calibration curves)
    solved_eq_c <- solving_equations(rv$fileimport_calibration,
                                     rv$choices_list,
                                     type = 1,
                                     rv = rv,
                                     mode = "corrected",
                                     logfilename = logfilename,
                                     minmax = rv$minmax)
    rv$fileimport_cal_corrected_c <- solved_eq_c[["results"]]
    colnames(rv$fileimport_cal_corrected_c) <- colnames(
      rv$fileimport_calibration
    )
    rv$substitutions_corrected_c <- solved_eq_c[["substitutions"]]

    expect_snapshot(
      x = round(table_prep(rv$fileimport_cal_corrected_c), 2),
      cran = FALSE,
      error = FALSE
    )
    expect_snapshot(
      x = round(table_prep(rv$substitutions_corrected_c), 2),
      cran = FALSE,
      error = FALSE
    )

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      rv$fileimport_cal_corrected_c,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      logfilename = logfilename,
      minmax = rv$minmax,
      seed = rv$seed
    )
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list_cubic <- testhelper_apply_robust_results_list(
      regression_results[["result_list"]]
    )
    # save regression statistics to reactive value
    rv$reg_stats_corrected_c <- statistics_list(rv$result_list_cubic,
                                                minmax = rv$minmax)

    expect_snapshot(
      x = round(table_prep(rv$reg_stats_corrected_c), 2),
      cran = FALSE,
      error = FALSE
    )

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

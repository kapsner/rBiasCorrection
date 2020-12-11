context("test functioning of algorithm, type 1")

prefix <- tempdir()
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)

test_that(
  desc = "algorithm test, type 1, minmax = TRUE selection_method = RelError",
  code = {

    local_edition(3)
    suppressWarnings(future::plan("multiprocess"))

    #"skip_on_cran()

    rv$minmax <- TRUE
    rv$selection_method <- "RelError"
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

    # reconstruct parts from app_plottingUtility.R
    regression_results <- regression_utility(rv$fileimport_calibration,
                                             "Testlocus",
                                             locus_id = NULL,
                                             rv = rv,
                                             mode = NULL,
                                             logfilename,
                                             minmax = rv$minmax,
                                             seed = rv$seed)
    plotlist_reg <- regression_results[["plot_list"]]
    rv$result_list <- regression_results[["result_list"]]

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

    # hyperbolic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = F
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
    rv$result_list_hyperbolic <- regression_results[["result_list"]]
    # save regression statistics to reactive value
    rv$reg_stats_corrected_h <- statistics_list(rv$result_list_hyperbolic,
                                                minmax = rv$minmax)

    # cubic correction
    rv$choices_list <- rv$reg_stats[, c("Name"), with = F
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
    rv$result_list_cubic <- regression_results[["result_list"]]
    # save regression statistics to reactive value
    rv$reg_stats_corrected_c <- statistics_list(rv$result_list_cubic,
                                                minmax = rv$minmax)

    # calculate final results
    rv$choices_list <- better_model(
      statstable_pre = rv$reg_stats,
      statstable_post_hyperbolic = rv$reg_stats_corrected_h,
      statstable_post_cubic = rv$reg_stats_corrected_c,
      selection_method = rv$selection_method
    )
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
    colnames(rv$fileimport_cal_corrected) <- colnames(
      rv$fileimport_calibration
    )

    # some tests
    expect_type(solved_eq, "list")
    expect_snapshot_value(
      x = solved_eq,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq, "cd7926e6e0")
    #   expect_known_hash(solved_eq, "0990d0bcd9")
    # }, class = "error", regexp = "cd7926e6e0|0990d0bcd9") # 0990d0bcd9,
    # # 13cff31610
    # nolint end
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_snapshot_value(
      x = rv$final_results,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$final_results, "47c0df77e1")
    #   expect_known_hash(rv$final_results, "93881d6d42")
    # }, class = "error", regexp = "47c0df77e1|93881d6d42") # 93881d6d42,
    # # eea1c59606
    # nolint end
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_snapshot_value(
      x = rv$substitutions,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$substitutions, "9b5110fcfe")
    #   expect_known_hash(rv$substitutions, "510026d492")
    # }, class = "error", regexp = "9b5110fcfe|510026d492") # 510026d492,
    # # 98e25743ed
    # nolint end
    expect_type(solved_eq2, "list")
    expect_snapshot_value(
      x = solved_eq2,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq2, "cf96c67f4c")
    #   expect_known_hash(solved_eq2, "a670800b3a")
    # }, class = "error", regexp = "cf96c67f4c|a670800b3a") # a670800b3a,
    # # 569f636794
    # nolint end
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_snapshot_value(
      x = rv$fileimport_cal_corrected,
      style = "serialize",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$fileimport_cal_corrected, "341a01cdf9")
    #   expect_known_hash(rv$fileimport_cal_corrected, "0347a748fe")
    # }, class = "error", regexp = "341a01cdf9|0347a748fe") # 0347a748fe,
    # # 5b8a8f6887
    # nolint end

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

context("test functioning of algorithm, type 1")

prefix <- tempdir()
# prefix <- "tests/testthat/" # nolint

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that(
  desc = "algorithm test, type 1, minmax = FALSE",
  code = {

    local_edition(3)

    suppressWarnings(future::plan("multiprocess"))

    #"skip_on_cran()

    rv$minmax <- FALSE
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

    # some tests
    expect_type(regression_results, "list")
    #" expect_known_hash(regression_results, "a75be8d5af")
    # oder 0bdeacf677, fc7ae30d08
    expect_type(plotlist_reg, "list")
    #" expect_known_hash(plotlist_reg, "20fa85b532")
    # oder c2e96f84fc, 0c3c5db52b

    expect_snapshot_value(
      x = rv$result_list,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$result_list, "66dbae3141")
    #   expect_known_hash(rv$result_list, "70c4aa7531")
    # }, class = "error", regexp = "66dbae3141|70c4aa7531")
    # # 70c4aa7531, d7f426a1a8, 2eb93a74d3
    # nolint end

    expect_type(rv$reg_stats, "list")
    expect_s3_class(rv$reg_stats, "data.table")

    expect_snapshot_value(
      x = rv$reg_stats,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$reg_stats, "3603e3abcb")
    #   expect_known_hash(rv$reg_stats, "d2e883e81d")
    # }, class = "error", regexp = "3603e3abcb|d2e883e81d") # d2e883e81d
    # #a27d84167e, b88f6a9fcf, 057c7d0a13, 33e1d855be
    # nolint end

    expect_equal(
      object = regression_results[["result_list"]],
      expected = regression_results2[["result_list"]],
      ignore_function_env = TRUE
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
    expect_type(solved_eq, "list")
    expect_snapshot_value(
      x = solved_eq,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq, "af9980684d")
    #   expect_known_hash(solved_eq, "15eb3be08d")
    # }, class = "error", regexp = "af9980684d|15eb3be08d") # 15eb3be08d
    # nolint end
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_snapshot_value(
      x = rv$final_results,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$final_results, "16edb2d3ca")
    #   expect_known_hash(rv$final_results, "c28704311f")
    # }, class = "error", regexp = "16edb2d3ca|c28704311f") # c28704311f
    # nolint end
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_snapshot_value(
      x = rv$substitutions,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$substitutions, "c67561996f")
    #   expect_known_hash(rv$substitutions, "6ca563f935")
    # }, class = "error", regexp = "c67561996f|6ca563f935") # 6ca563f935
    # nolint end
    expect_type(solved_eq2, "list")
    expect_snapshot_value(
      x = solved_eq2,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq2, "9ea9e3d04d")
    #   expect_known_hash(solved_eq2, "b4605bb70f")
    # }, class = "error", regexp = "9ea9e3d04d|b4605bb70f") # b4605bb70f
    # nolint end
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_snapshot_value(
      x = rv$fileimport_cal_corrected,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$fileimport_cal_corrected, "947a883caa")
    #   expect_known_hash(rv$fileimport_cal_corrected, "fd7708da93")
    # }, class = "error", regexp = "947a883caa|fd7708da93") # fd7708da93
    # nolint end



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

    expect_type(solved_eq_h, "list")
    expect_snapshot_value(
      x = solved_eq_h,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq_h, "8ce2d1e597")
    #   expect_known_hash(solved_eq_h, "522a653888")
    # }, class = "error", regexp = "8ce2d1e597|522a653888") # 522a653888
    # nolint end
    expect_type(rv$fileimport_cal_corrected_h, "list")
    expect_s3_class(rv$fileimport_cal_corrected_h, "data.table")
    expect_snapshot_value(
      x = rv$fileimport_cal_corrected_h,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$fileimport_cal_corrected_h, "0592bdf33f")
    #   expect_known_hash(rv$fileimport_cal_corrected_h, "652ac3a81f")
    # }, class = "error", regexp = "0592bdf33f|652ac3a81f") # 652ac3a81f
    # nolint end
    expect_type(rv$substitutions_corrected_h, "list")
    expect_s3_class(rv$substitutions_corrected_h, "data.table")
    expect_snapshot_value(
      x = rv$substitutions_corrected_h,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$substitutions_corrected_h, "d1b7a8b151")
    #   expect_known_hash(rv$substitutions_corrected_h, "16c7efb563")
    # }, class = "error", regexp = "d1b7a8b151|16c7efb563") # 16c7efb563
    # nolint end

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

    expect_type(regression_results, "list")
    #" expect_known_hash(regression_results, "a75be8d5af")
    # oder 0bdeacf677, fc7ae30d08
    expect_type(plotlist_reg, "list")
    #" expect_known_hash(plotlist_reg, "20fa85b532")
    # oder c2e96f84fc, 0c3c5db52b
    expect_type(rv$result_list_hyperbolic, "list")
    expect_snapshot_value(
      x = rv$result_list_hyperbolic,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$result_list_hyperbolic, "b39e4cbca1")
    #   expect_known_hash(rv$result_list_hyperbolic, "f49cd1b17c")
    # }, class = "error", regexp = "b39e4cbca1|f49cd1b17c") # f49cd1b17c
    # 3d50611917, ccbe9ff93a
    # nolint end
    expect_type(rv$reg_stats_corrected_h, "list")
    expect_s3_class(rv$reg_stats_corrected_h, "data.table")
    expect_snapshot_value(
      x = rv$reg_stats_corrected_h,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$reg_stats_corrected_h, "f8ac0ddffd")
    #   expect_known_hash(rv$reg_stats_corrected_h, "9b4f5b4dcb")
    # }, class = "error", regexp = "f8ac0ddffd|9b4f5b4dcb") # 9b4f5b4dcb
    # #a27d84167e, 5205aae446, 81e6bcb79b, 5c8f64e551
    # nolint end


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

    expect_type(solved_eq_c, "list")
    expect_snapshot_value(
      x = solved_eq_c,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq_c, "a840531423")
    #   expect_known_hash(solved_eq_c, "b708e60ff0")
    # }, class = "error", regexp = "a840531423|b708e60ff0") # b708e60ff0
    # nolint end
    expect_type(rv$fileimport_cal_corrected_c, "list")
    expect_s3_class(rv$fileimport_cal_corrected_c, "data.table")
    expect_snapshot_value(
      x = rv$fileimport_cal_corrected_c,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$fileimport_cal_corrected_c, "a99f550089")
    #   expect_known_hash(rv$fileimport_cal_corrected_c, "9c02bbabfd")
    # }, class = "error", regexp = "a99f550089|9c02bbabfd") # 9c02bbabfd
    # nolint end
    expect_type(rv$substitutions_corrected_c, "list")
    expect_s3_class(rv$substitutions_corrected_c, "data.table")
    expect_snapshot_value(
      x = rv$substitutions_corrected_c,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_known_hash(rv$substitutions_corrected_c, "5e15c67e45")
    # nolint end

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

    expect_type(regression_results, "list")
    #" expect_known_hash(regression_results, "a75be8d5af")
    # oder 0bdeacf677, fc7ae30d08
    expect_type(plotlist_reg, "list")
    #" expect_known_hash(plotlist_reg, "20fa85b532")
    # oder c2e96f84fc, 0c3c5db52b
    expect_type(rv$result_list_cubic, "list")
    expect_snapshot_value(
      x = rv$result_list_cubic,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$result_list_cubic, "77dda04aba")
    #   expect_known_hash(rv$result_list_cubic, "36e27655ea")
    # }, class = "error", regexp = "77dda04aba|36e27655ea") # 36e27655ea
    # # 7214d93552, 9bc037ad08
    # nolint end

    expect_type(rv$reg_stats_corrected_c, "list")
    expect_s3_class(rv$reg_stats_corrected_c, "data.table")
    expect_snapshot_value(
      x = rv$reg_stats_corrected_c,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$reg_stats_corrected_c, "446af74c7a")
    #   expect_known_hash(rv$reg_stats_corrected_c, "88bb1a0300")
    # }, class = "error", regexp = "446af74c7a|88bb1a0300") # 88bb1a0300
    # #a27d84167e, 1d48c373f6, 90a3a2cb09, e79434dab7
    # nolint end
  })


test_that(
  desc = paste0("algorithm test, type 1, minmax = FALSE ",
                "selection_method = RelError"),
  code = {

    local_edition(3)
    suppressWarnings(future::plan("multiprocess"))

    #"skip_on_cran()

    rv$minmax <- FALSE
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

    regression_results2 <- regression_type1(
      rv$fileimport_calibration,
      rv$vec_cal,
      mode = NULL,
      logfilename,
      minmax = rv$minmax,
      locus_id = NULL,
      locusname = rv$sample_locus_name,
      seed = rv$seed
    )

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

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      rv$fileimport_cal_corrected_c,
      samplelocusname = rv$sample_locus_name,
      rv = rv, mode = "corrected",
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
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq, "d9350cb0b2")
    #   expect_known_hash(solved_eq, "fe0adc367e")
    # }, class = "error", regexp = "d9350cb0b2|fe0adc367e") # fe0adc367e
    # nolint end
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_snapshot_value(
      x = rv$final_results,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$final_results, "3189653abe")
    #   expect_known_hash(rv$final_results, "b25fab0cf3")
    # }, class = "error", regexp = "3189653abe|b25fab0cf3") # b25fab0cf3
    # nolint end
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_snapshot_value(
      x = rv$substitutions,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$substitutions, "ca18ba24a0")
    #   expect_known_hash(rv$substitutions, "7e1aaead80")
    # }, class = "error", regexp = "ca18ba24a0|7e1aaead80") # 7e1aaead80
    # nolint end
    expect_type(solved_eq2, "list")
    expect_snapshot_value(
      x = solved_eq2,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(solved_eq2, "3e2fa23575")
    #   expect_known_hash(solved_eq2, "df59cd9f96")
    # }, class = "error", regexp = "3e2fa23575|df59cd9f96") # df59cd9f96
    # nolint end
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_snapshot_value(
      x = rv$fileimport_cal_corrected,
      style = "serialize",
      cran = FALSE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(rv$fileimport_cal_corrected, "ff0a502640")
    #   expect_known_hash(rv$fileimport_cal_corrected, "248f22f000")
    # }, class = "error", regexp = "ff0a502640|248f22f000") # 248f22f000
    # nolint end

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

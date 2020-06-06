context("test functioning of algorithm, type 1")

prefix <- "./"
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that(
  desc = "algorithm test, type 1, minmax = FALSE",
  code = {

    #"skip_on_cran()

    rv$minmax <- FALSE
    rv$sample_locus_name <- "Test"
    rv$seed <- 1234

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

    # some tests
    expect_length(rv$vec_cal, 10)
    expect_type(rv$vec_cal, "character")


    # reconstruct parts from app_plottingUtility.R
    regression_results <- regression_utility(rv$fileimport_calibration,
                                             "Testlocus",
                                             locus_id = NULL,
                                             rv = rv,
                                             mode = NULL,
                                             headless = TRUE,
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

    # some tests
    expect_type(regression_results, "list")
    #" expect_known_hash(regression_results, "a75be8d5af")
    # oder 0bdeacf677, fc7ae30d08
    expect_type(plotlist_reg, "list")
    #" expect_known_hash(plotlist_reg, "20fa85b532")
    # oder c2e96f84fc, 0c3c5db52b
    #"expect_known_hash(rv$result_list, "66dbae3141")
    # 70c4aa7531, d7f426a1a8, 2eb93a74d3
    expect_known_hash(unlist(rv$result_list[["row_means"]]), "645f15f793")

    expect_type(rv$reg_stats, "list")
    expect_s3_class(rv$reg_stats, "data.table")
    expect_known_hash(rv$reg_stats, "3603e3abcb") # d2e883e81d
    #a27d84167e, b88f6a9fcf, 057c7d0a13, 33e1d855be
    expect_equal(regression_results, regression_results2)
    expect_equal(regression_results[["plot_list"]],
                 regression_results2[["plot_list"]])
    expect_equal(regression_results[["result_list"]],
                 regression_results2[["result_list"]])

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
    expect_known_hash(solved_eq, "af9980684d") # 15eb3be08d
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_known_hash(rv$final_results, "16edb2d3ca") # c28704311f
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_known_hash(rv$substitutions, "c67561996f") # 6ca563f935
    expect_type(solved_eq2, "list")
    expect_known_hash(solved_eq2, "9ea9e3d04d") # b4605bb70f
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected, "947a883caa")
    # fd7708da93



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
    expect_known_hash(solved_eq_h, "8ce2d1e597") # 522a653888
    expect_type(rv$fileimport_cal_corrected_h, "list")
    expect_s3_class(rv$fileimport_cal_corrected_h, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected_h, "0592bdf33f")
    # 652ac3a81f
    expect_type(rv$substitutions_corrected_h, "list")
    expect_s3_class(rv$substitutions_corrected_h, "data.table")
    expect_known_hash(rv$substitutions_corrected_h, "d1b7a8b151")
    # 16c7efb563

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      rv$fileimport_cal_corrected_h,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      headless = TRUE,
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
    #"expect_known_hash(rv$result_list_hyperbolic, "b39e4cbca1") # f49cd1b17c
    # 3d50611917, ccbe9ff93a
    expect_known_hash(
      unlist(rv$result_list_hyperbolic[["row_means"]]),
      "5144e200a9"
    )
    expect_type(rv$reg_stats_corrected_h, "list")
    expect_s3_class(rv$reg_stats_corrected_h, "data.table")
    #"expect_known_hash(rv$reg_stats_corrected_h, "f8ac0ddffd") # 9b4f5b4dcb
    #a27d84167e, 5205aae446, 81e6bcb79b, 5c8f64e551
    expect_known_hash(
      as.matrix(rv$reg_stats_corrected_h[10, ]),
      "bc0545975d"
    )


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
    expect_known_hash(solved_eq_c, "a840531423") # b708e60ff0
    expect_type(rv$fileimport_cal_corrected_c, "list")
    expect_s3_class(rv$fileimport_cal_corrected_c, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected_c, "a99f550089")
    # 9c02bbabfd
    expect_type(rv$substitutions_corrected_c, "list")
    expect_s3_class(rv$substitutions_corrected_c, "data.table")
    expect_known_hash(rv$substitutions_corrected_c, "5e15c67e45")

    # calculate new calibration curves from corrected calibration data
    regression_results <- regression_utility(
      rv$fileimport_cal_corrected_c,
      samplelocusname = rv$sample_locus_name,
      rv = rv,
      mode = "corrected",
      headless = TRUE,
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
    #"expect_known_hash(rv$result_list_cubic, "77dda04aba") # 36e27655ea
    # 7214d93552, 9bc037ad08
    expect_known_hash(
      unlist(rv$result_list_cubic[["row_means"]]),
      "c39c72b05c"
    )

    expect_type(rv$reg_stats_corrected_c, "list")
    expect_s3_class(rv$reg_stats_corrected_c, "data.table")
    #"expect_known_hash(rv$reg_stats_corrected_c, "446af74c7a") # 88bb1a0300
    #a27d84167e, 1d48c373f6, 90a3a2cb09, e79434dab7
    expect_known_hash(
      as.matrix(rv$reg_stats_corrected_c[10, ]),
      "6900d1046f"
    )
  })


test_that(
  desc = paste0("algorithm test, type 1, minmax = FALSE ",
                "selection_method = RelError"),
  code = {

    #"skip_on_cran()

    rv$minmax <- FALSE
    rv$selection_method <- "RelError"
    rv$sample_locus_name <- "Test"
    rv$seed <- 1234

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
    regression_results <- regression_utility(rv$fileimport_calibration,
                                             "Testlocus",
                                             locus_id = NULL,
                                             rv = rv,
                                             mode = NULL,
                                             headless = TRUE,
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
      headless = TRUE,
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
      headless = TRUE,
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
    expect_known_hash(solved_eq, "d9350cb0b2") # fe0adc367e
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_known_hash(rv$final_results, "3189653abe") # b25fab0cf3
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_known_hash(rv$substitutions, "ca18ba24a0") # 7e1aaead80
    expect_type(solved_eq2, "list")
    expect_known_hash(solved_eq2, "3e2fa23575") # df59cd9f96
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected, "ff0a502640")
    # 248f22f000
  })

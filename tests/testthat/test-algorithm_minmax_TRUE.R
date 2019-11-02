context("test functioning of algorithm, type 1")

prefix <- "./"
#" prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that(
  desc = "algorithm test, type 1, minmax = TRUE",
  code = {

    skip_on_cran()

    rv$minmax <- TRUE
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
    expect_type(rv$result_list, "list")
    expect_known_hash(rv$result_list, "e23adb8bba")
    # 4c700dcb89, 8c7d29964f
    expect_type(rv$reg_stats, "list")
    expect_s3_class(rv$reg_stats, "data.table")
    expect_known_hash(rv$reg_stats, "79c54720c8") # 261933672d
    #f20afd797f, 8d9e9c577f, 3b8a2c9335, a416f49f04
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
    expect_known_hash(solved_eq, "9ba9164e21") # e17032de5b, eb472ae0e7
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_known_hash(rv$final_results, "b4e0e01ac2") # 16fca713ba, 920658389f
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_known_hash(rv$substitutions, "161577b615") # 6f50a58a2f, 30c692c633
    expect_type(solved_eq2, "list")
    expect_known_hash(solved_eq2, "ef6f39b1ce") # 202891fe22, dd15288aba
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected, "215eb69643")
    # d082d296f8, 913f716d0c



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
    expect_known_hash(solved_eq_h, "2a807cdbb3") # 8a58318532, ac4aee295b
    expect_type(rv$fileimport_cal_corrected_h, "list")
    expect_s3_class(rv$fileimport_cal_corrected_h, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected_h, "8d01172516")
    # 5e110ecf0d, fc5617597f
    expect_type(rv$substitutions_corrected_h, "list")
    expect_s3_class(rv$substitutions_corrected_h, "data.table")
    expect_known_hash(rv$substitutions_corrected_h, "33afa269a4")

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
    expect_known_hash(rv$result_list_hyperbolic, "9c70512014")
    # 0f3e987b00, 52ce26f8c5
    expect_type(rv$reg_stats_corrected_h, "list")
    expect_s3_class(rv$reg_stats_corrected_h, "data.table")
    expect_known_hash(rv$reg_stats_corrected_h, "f4b8df4ad1") # 87bf0a0b86
    #e128ff333d, aa7217b008, 22990dacfc, 46ca43f245



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
    expect_known_hash(solved_eq_c, "fa834b6f83") # dbac3589ca, a840531423
    expect_type(rv$fileimport_cal_corrected_c, "list")
    expect_s3_class(rv$fileimport_cal_corrected_c, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected_c, "a42c8c2f2d")
    # d0447e2521, a99f550089
    expect_type(rv$substitutions_corrected_c, "list")
    expect_s3_class(rv$substitutions_corrected_c, "data.table")
    expect_known_hash(rv$substitutions_corrected_c, "3e2bca3b0a")
    # 33afa269a4, 5e15c67e45

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
    expect_known_hash(rv$result_list_cubic, "209f8c844d")
    # dcd8ba3827, 7c9569a4a4
    expect_type(rv$reg_stats_corrected_c, "list")
    expect_s3_class(rv$reg_stats_corrected_c, "data.table")
    expect_known_hash(rv$reg_stats_corrected_c, "09ad550c5e") # c38ea3ed70
    #b41b6cc539, fe5ea3da9c, 050face677, 4fb40e13ea
  })

test_that(
  desc = "algorithm test, type 1, minmax = TRUE selection_method = RelError",
  code = {

    skip_on_cran()

    rv$minmax <- TRUE
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
    rv$substitutions_corrected_c <- solved_eq_c[["substitutions"]]

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
    expect_known_hash(solved_eq, "0990d0bcd9") # cd7926e6e0, 13cff31610
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_known_hash(rv$final_results, "93881d6d42") # 47c0df77e1, eea1c59606
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_known_hash(rv$substitutions, "510026d492") # 9b5110fcfe, 98e25743ed
    expect_type(solved_eq2, "list")
    expect_known_hash(solved_eq2, "a670800b3a") # cf96c67f4c, 569f636794
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_known_hash(rv$fileimport_cal_corrected, "0347a748fe")
    # 341a01cdf9, 5b8a8f6887
  })

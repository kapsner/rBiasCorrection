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

    #"skip_on_cran()

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
    expect_error({
      expect_known_hash(rv$result_list, "4c700dcb89")
      expect_known_hash(rv$result_list, "e23adb8bba")
    }, class = "error", regexp = "4c700dcb89|e23adb8bba") # e23adb8bba
    # 8c7d29964f

    expect_type(rv$reg_stats, "list")
    expect_s3_class(rv$reg_stats, "data.table")
    expect_error({
      expect_known_hash(rv$reg_stats, "261933672d")
      expect_known_hash(rv$reg_stats, "79c54720c8")
    }, class = "error", regexp = "261933672d|79c54720c8") # 79c54720c8
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
    expect_error({
      expect_known_hash(solved_eq, "0469f8831e")
      expect_known_hash(solved_eq, "9ba9164e21")
    }, class = "error", regexp = "0469f8831e|9ba9164e21") # e17032de5b
    # 9ba9164e21, eb472ae0e7
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_error({
      expect_known_hash(rv$final_results, "fb751fd42e")
      expect_known_hash(rv$final_results, "b4e0e01ac2")
    }, class = "error", regexp = "fb751fd42e|b4e0e01ac2") # 16fca713ba
    # b4e0e01ac2, 920658389f
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_error({
      expect_known_hash(rv$substitutions, "6f50a58a2f")
      expect_known_hash(rv$substitutions, "161577b615")
    }, class = "error", regexp = "6f50a58a2f|161577b615") # 161577b615,
    # 30c692c633
    expect_type(solved_eq2, "list")
    expect_error({
      expect_known_hash(solved_eq2, "202891fe22")
      expect_known_hash(solved_eq2, "ef6f39b1ce")
    }, class = "error", regexp = "202891fe22|ef6f39b1ce") # ef6f39b1ce,
    # dd15288aba
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_error({
      expect_known_hash(rv$fileimport_cal_corrected, "d082d296f8")
      expect_known_hash(rv$fileimport_cal_corrected, "215eb69643")
    }, class = "error", regexp = "d082d296f8|215eb69643") # 215eb69643,
    # 913f716d0c



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
    expect_error({
      expect_known_hash(solved_eq_h, "8a58318532")
      expect_known_hash(solved_eq_h, "2a807cdbb3")
    }, class = "error", regexp = "8a58318532|2a807cdbb3") # 2a807cdbb3,
    # ac4aee295b
    expect_type(rv$fileimport_cal_corrected_h, "list")
    expect_s3_class(rv$fileimport_cal_corrected_h, "data.table")
    expect_error({
      expect_known_hash(rv$fileimport_cal_corrected_h, "5e110ecf0d")
      expect_known_hash(rv$fileimport_cal_corrected_h, "8d01172516")
    }, class = "error", regexp = "5e110ecf0d|8d01172516") # 8d01172516,
    # fc5617597f
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
    expect_error({
      expect_known_hash(rv$result_list_hyperbolic, "0f3e987b00")
      expect_known_hash(rv$result_list_hyperbolic, "9c70512014")
    }, class = "error", regexp = "0f3e987b00|9c70512014") # 9c70512014,
    # 52ce26f8c5

    expect_type(rv$reg_stats_corrected_h, "list")
    expect_s3_class(rv$reg_stats_corrected_h, "data.table")
    expect_error({
      expect_known_hash(rv$reg_stats_corrected_h, "87bf0a0b86")
      expect_known_hash(rv$reg_stats_corrected_h, "f4b8df4ad1")
    }, class = "error", regexp = "87bf0a0b86|f4b8df4ad1") # f4b8df4ad1
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
    expect_error({
      expect_known_hash(solved_eq_c, "dbac3589ca")
      expect_known_hash(solved_eq_c, "fa834b6f83")
    }, class = "error", regexp = "dbac3589ca|fa834b6f83") # fa834b6f83,
    # a840531423
    expect_type(rv$fileimport_cal_corrected_c, "list")
    expect_s3_class(rv$fileimport_cal_corrected_c, "data.table")
    expect_error({
      expect_known_hash(rv$fileimport_cal_corrected_c, "d0447e2521")
      expect_known_hash(rv$fileimport_cal_corrected_c, "a42c8c2f2d")
    }, class = "error", regexp = "d0447e2521|a42c8c2f2d") # a42c8c2f2d,
    # a99f550089
    expect_type(rv$substitutions_corrected_c, "list")
    expect_s3_class(rv$substitutions_corrected_c, "data.table")
    expect_error({
      expect_known_hash(rv$substitutions_corrected_c, "33afa269a4")
      expect_known_hash(rv$substitutions_corrected_c, "3e2bca3b0a")
    }, class = "error", regexp = "33afa269a4|3e2bca3b0a") # 3e2bca3b0a,
    # 5e15c67e45

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
    expect_error({
      expect_known_hash(rv$result_list_cubic, "dcd8ba3827")
      expect_known_hash(rv$result_list_cubic, "209f8c844d")
    }, class = "error", regexp = "dcd8ba3827|209f8c844d") # 209f8c844d,
    # 7c9569a4a4

    expect_type(rv$reg_stats_corrected_c, "list")
    expect_s3_class(rv$reg_stats_corrected_c, "data.table")
    expect_error({
      expect_known_hash(rv$reg_stats_corrected_c, "c38ea3ed70")
      expect_known_hash(rv$reg_stats_corrected_c, "09ad550c5e")
    }, class = "error", regexp = "c38ea3ed70|09ad550c5e") # 09ad550c5e
    #b41b6cc539, fe5ea3da9c, 050face677, 4fb40e13ea
  })

test_that(
  desc = "algorithm test, type 1, minmax = TRUE selection_method = RelError",
  code = {

    #"skip_on_cran()

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
    expect_error({
      expect_known_hash(solved_eq, "cd7926e6e0")
      expect_known_hash(solved_eq, "0990d0bcd9")
    }, class = "error", regexp = "cd7926e6e0|0990d0bcd9") # 0990d0bcd9,
    # 13cff31610
    expect_type(rv$final_results, "list")
    expect_s3_class(rv$final_results, "data.table")
    expect_error({
      expect_known_hash(rv$final_results, "47c0df77e1")
      expect_known_hash(rv$final_results, "93881d6d42")
    }, class = "error", regexp = "47c0df77e1|93881d6d42") # 93881d6d42,
    # eea1c59606
    expect_type(rv$substitutions, "list")
    expect_s3_class(rv$substitutions, "data.table")
    expect_error({
      expect_known_hash(rv$substitutions, "9b5110fcfe")
      expect_known_hash(rv$substitutions, "510026d492")
    }, class = "error", regexp = "9b5110fcfe|510026d492") # 510026d492,
    # 98e25743ed
    expect_type(solved_eq2, "list")
    expect_error({
      expect_known_hash(solved_eq2, "cf96c67f4c")
      expect_known_hash(solved_eq2, "a670800b3a")
    }, class = "error", regexp = "cf96c67f4c|a670800b3a") # a670800b3a,
    # 569f636794
    expect_type(rv$fileimport_cal_corrected, "list")
    expect_s3_class(rv$fileimport_cal_corrected, "data.table")
    expect_error({
      expect_known_hash(rv$fileimport_cal_corrected, "341a01cdf9")
      expect_known_hash(rv$fileimport_cal_corrected, "0347a748fe")
    }, class = "error", regexp = "341a01cdf9|0347a748fe") # 0347a748fe,
    # 5b8a8f6887
  })

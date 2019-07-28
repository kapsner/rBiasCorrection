context("test functioning of algorithm, type 1")

prefix <- "./"
#prefix <- "tests/testthat/"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")

# initialize our list for reactive values
rv <- list()


library(data.table)


test_that("algorithm test, type 1, minmax = TRUE",{
  rv$minmax <- TRUE

  # experimental data
  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  rv$fileimportExp <- cleanDT_(exp_type_1, "experimental", 1, logfilename)[["dat"]]

  # calibration data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
  cal_type_1 <- cleanDT_(cal_type_1, "calibration", 1, logfilename)
  rv$fileimportCal <- cal_type_1[["dat"]]
  rv$vec_cal <- cal_type_1[["vec_cal"]]

  # some tests
  expect_length(rv$vec_cal, 10)
  expect_type(rv$vec_cal, "character")


  # reconstruct parts from app_plottingUtility.R
  regression_results <- regressionUtility_(rv$fileimportCal, "Testlocus", locus_id = NULL, rv = rv, mode = NULL, headless = TRUE, logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list <- regression_results[["result_list"]]

  regression_results2 <- regression_type1(rv$fileimportCal, rv$vec_cal, mode=NULL, logfilename, minmax = rv$minmax)

  # save regression statistics to reactive value
  rv$regStats <- statisticsList_(rv$result_list, minmax = rv$minmax)

  # some tests
  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list, "list")
  expect_known_hash(rv$result_list, "8c7d29964f")
  expect_type(rv$regStats, "list")
  expect_s3_class(rv$regStats, "data.table")
  expect_known_hash(rv$regStats, "8d9e9c577f") #f20afd797f
  expect_equal(regression_results, regression_results2)
  expect_equal(regression_results[["plot_list"]], regression_results2[["plot_list"]])
  expect_equal(regression_results[["result_list"]], regression_results2[["result_list"]])

  # calculate final results
  # default rv$choices_list == rv$regStats[,.(Name, better_model)]
  solved_eq <- solvingEquations_(rv$fileimportExp, rv$regStats[,c("Name", "better_model"),with=F], type = 1, rv = rv, logfilename = logfilename, minmax = rv$minmax)
  rv$finalResults <- solved_eq[["results"]]
  rv$substitutions <- solved_eq[["substitutions"]]

  # Calibration Data (to show corrected calibration curves)
  solved_eq2 <- solvingEquations_(rv$fileimportCal, rv$regStats[,c("Name", "better_model"),with=F], type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected <- solved_eq2[["results"]]
  colnames(rv$fileimportCal_corrected) <- colnames(rv$fileimportCal)

  # some tests
  expect_type(solved_eq, "list")
  expect_known_hash(solved_eq, "564019ef97")
  expect_type(rv$finalResults, "list")
  expect_s3_class(rv$finalResults, "data.table")
  expect_known_hash(rv$finalResults, "920658389f")
  expect_type(rv$substitutions, "list")
  expect_s3_class(rv$substitutions, "data.table")
  expect_known_hash(rv$substitutions, "411246447a")
  expect_type(solved_eq2, "list")
  expect_known_hash(solved_eq2, "91795b9115")
  expect_type(rv$fileimportCal_corrected, "list")
  expect_s3_class(rv$fileimportCal_corrected, "data.table")
  expect_known_hash(rv$fileimportCal_corrected, "913f716d0c")



  # hyperbolic correction
  rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=0]

  # correct calibration data (to show corrected calibration curves)
  solved_eq_h <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected_h <- solved_eq_h[["results"]]
  colnames(rv$fileimportCal_corrected_h) <- colnames(rv$fileimportCal)

  expect_type(solved_eq_h, "list")
  expect_known_hash(solved_eq_h, "70377bc9ca")
  expect_type(rv$fileimportCal_corrected_h, "list")
  expect_s3_class(rv$fileimportCal_corrected_h, "data.table")
  expect_known_hash(rv$fileimportCal_corrected_h, "fc5617597f")

  # calculate new calibration curves from corrected calibration data
  regression_results <- regressionUtility_(rv$fileimportCal_corrected_h, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list_hyperbolic <- regression_results[["result_list"]]
  # save regression statistics to reactive value
  rv$regStats_corrected_h <- statisticsList_(rv$result_list_hyperbolic, minmax = rv$minmax)

  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list_hyperbolic, "list")
  expect_known_hash(rv$result_list_hyperbolic, "52ce26f8c5")
  expect_type(rv$regStats_corrected_h, "list")
  expect_s3_class(rv$regStats_corrected_h, "data.table")
  expect_known_hash(rv$regStats_corrected_h, "aa7217b008") #e128ff333d



  # cubic correction
  rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=1]

  # correct calibration data (to show corrected calibration curves)
  solved_eq_c <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected_c <- solved_eq_c[["results"]]
  colnames(rv$fileimportCal_corrected_c) <- colnames(rv$fileimportCal)

  expect_type(solved_eq_c, "list")
  expect_known_hash(solved_eq_c, "4bced23098")
  expect_type(rv$fileimportCal_corrected_c, "list")
  expect_s3_class(rv$fileimportCal_corrected_c, "data.table")
  expect_known_hash(rv$fileimportCal_corrected_c, "a99f550089")

  # calculate new calibration curves from corrected calibration data
  regression_results <- regressionUtility_(rv$fileimportCal_corrected_c, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list_cubic <- regression_results[["result_list"]]
  # save regression statistics to reactive value
  rv$regStats_corrected_c <- statisticsList_(rv$result_list_cubic, minmax = rv$minmax)

  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list_cubic, "list")
  expect_known_hash(rv$result_list_cubic, "7c9569a4a4")
  expect_type(rv$regStats_corrected_c, "list")
  expect_s3_class(rv$regStats_corrected_c, "data.table")
  expect_known_hash(rv$regStats_corrected_c, "fe5ea3da9c") #b41b6cc539
})


test_that("algorithm test, type 1, minmax = FALSE",{
  rv$minmax <- FALSE

  # experimental data
  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  rv$fileimportExp <- cleanDT_(exp_type_1, "experimental", 1, logfilename)[["dat"]]

  # calibration data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
  cal_type_1 <- cleanDT_(cal_type_1, "calibration", 1, logfilename)
  rv$fileimportCal <- cal_type_1[["dat"]]
  rv$vec_cal <- cal_type_1[["vec_cal"]]

  # some tests
  expect_length(rv$vec_cal, 10)
  expect_type(rv$vec_cal, "character")


  # reconstruct parts from app_plottingUtility.R
  regression_results <- regressionUtility_(rv$fileimportCal, "Testlocus", locus_id = NULL, rv = rv, mode = NULL, headless = TRUE, logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list <- regression_results[["result_list"]]

  regression_results2 <- regression_type1(rv$fileimportCal, rv$vec_cal, mode=NULL, logfilename, minmax = rv$minmax)

  # save regression statistics to reactive value
  rv$regStats <- statisticsList_(rv$result_list, minmax = rv$minmax)

  # some tests
  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list, "list")
  expect_known_hash(rv$result_list, "d7f426a1a8")
  expect_type(rv$regStats, "list")
  expect_s3_class(rv$regStats, "data.table")
  expect_known_hash(rv$regStats, "b88f6a9fcf") #a27d84167e
  expect_equal(regression_results, regression_results2)
  expect_equal(regression_results[["plot_list"]], regression_results2[["plot_list"]])
  expect_equal(regression_results[["result_list"]], regression_results2[["result_list"]])

  # calculate final results
  # default rv$choices_list == rv$regStats[,.(Name, better_model)]
  solved_eq <- solvingEquations_(rv$fileimportExp, rv$regStats[,c("Name", "better_model"),with=F], type = 1, rv = rv, logfilename = logfilename, minmax = rv$minmax)
  rv$finalResults <- solved_eq[["results"]]
  rv$substitutions <- solved_eq[["substitutions"]]

  # Calibration Data (to show corrected calibration curves)
  solved_eq2 <- solvingEquations_(rv$fileimportCal, rv$regStats[,c("Name", "better_model"),with=F], type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected <- solved_eq2[["results"]]
  colnames(rv$fileimportCal_corrected) <- colnames(rv$fileimportCal)

  # some tests
  expect_type(solved_eq, "list")
  expect_known_hash(solved_eq, "f2112fc3f2")
  expect_type(rv$finalResults, "list")
  expect_s3_class(rv$finalResults, "data.table")
  expect_known_hash(rv$finalResults, "16edb2d3ca")
  expect_type(rv$substitutions, "list")
  expect_s3_class(rv$substitutions, "data.table")
  expect_known_hash(rv$substitutions, "910b14b78c")
  expect_type(solved_eq2, "list")
  expect_known_hash(solved_eq2, "d7540c0ec6")
  expect_type(rv$fileimportCal_corrected, "list")
  expect_s3_class(rv$fileimportCal_corrected, "data.table")
  expect_known_hash(rv$fileimportCal_corrected, "947a883caa")



  # hyperbolic correction
  rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=0]

  # correct calibration data (to show corrected calibration curves)
  solved_eq_h <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected_h <- solved_eq_h[["results"]]
  colnames(rv$fileimportCal_corrected_h) <- colnames(rv$fileimportCal)

  expect_type(solved_eq_h, "list")
  expect_known_hash(solved_eq_h, "7b7b33f9e1")
  expect_type(rv$fileimportCal_corrected_h, "list")
  expect_s3_class(rv$fileimportCal_corrected_h, "data.table")
  expect_known_hash(rv$fileimportCal_corrected_h, "0592bdf33f")

  # calculate new calibration curves from corrected calibration data
  regression_results <- regressionUtility_(rv$fileimportCal_corrected_h, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list_hyperbolic <- regression_results[["result_list"]]
  # save regression statistics to reactive value
  rv$regStats_corrected_h <- statisticsList_(rv$result_list_hyperbolic, minmax = rv$minmax)

  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list_hyperbolic, "list")
  expect_known_hash(rv$result_list_hyperbolic, "3d50611917")
  expect_type(rv$regStats_corrected_h, "list")
  expect_s3_class(rv$regStats_corrected_h, "data.table")
  expect_known_hash(rv$regStats_corrected_h, "5205aae446") #a27d84167e



  # cubic correction
  rv$choices_list <- rv$regStats[,c("Name"), with=F][,("better_model"):=1]

  # correct calibration data (to show corrected calibration curves)
  solved_eq_c <- solvingEquations_(rv$fileimportCal, rv$choices_list, type = 1, rv = rv, mode = "corrected", logfilename = logfilename, minmax = rv$minmax)
  rv$fileimportCal_corrected_c <- solved_eq_c[["results"]]
  colnames(rv$fileimportCal_corrected_c) <- colnames(rv$fileimportCal)

  expect_type(solved_eq_c, "list")
  expect_known_hash(solved_eq_c, "4bced23098")
  expect_type(rv$fileimportCal_corrected_c, "list")
  expect_s3_class(rv$fileimportCal_corrected_c, "data.table")
  expect_known_hash(rv$fileimportCal_corrected_c, "a99f550089")

  # calculate new calibration curves from corrected calibration data
  regression_results <- regressionUtility_(rv$fileimportCal_corrected_c, samplelocusname=rv$sampleLocusName, rv=rv, mode="corrected", headless = TRUE, logfilename = logfilename, minmax = rv$minmax)
  plotlistR <- regression_results[["plot_list"]]
  rv$result_list_cubic <- regression_results[["result_list"]]
  # save regression statistics to reactive value
  rv$regStats_corrected_c <- statisticsList_(rv$result_list_cubic, minmax = rv$minmax)

  expect_type(regression_results, "list")
  #expect_known_hash(regression_results, "a75be8d5af") # oder 0bdeacf677, fc7ae30d08
  expect_type(plotlistR, "list")
  #expect_known_hash(plotlistR, "20fa85b532") # oder c2e96f84fc, 0c3c5db52b
  expect_type(rv$result_list_cubic, "list")
  expect_known_hash(rv$result_list_cubic, "7214d93552")
  expect_type(rv$regStats_corrected_c, "list")
  expect_s3_class(rv$regStats_corrected_c, "data.table")
  expect_known_hash(rv$regStats_corrected_c, "1d48c373f6") #a27d84167e
})

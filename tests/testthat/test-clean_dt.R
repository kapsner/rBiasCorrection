context("dataimport filechecks")

prefix <- "./"

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "log.txt")

library(data.table)

test_that(
  desc = "test normal function of file import of type 1",
  code = {
  # experimental data
  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
  expect_known_hash(exp_type_1, "ab5287b084")

  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)
  expect_type(exp_type_1, "list")

  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1_empty_col.csv"),
                      header = T)
  exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
  expect_known_hash(exp_type_1, "ab5287b084")

  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1_empty_row.csv"),
                      header = T)
  exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
  expect_known_hash(exp_type_1, "1fb62d8498")

  # calibration data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
  cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
  expect_known_hash(cal_type_1, "23f21fc354")
})

test_that(
  desc = "test normal function of file import of type 2",
  code = {
  # experimental data
  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2.csv"))
  exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
  expect_known_hash(exp_type_2, "f8d57b6e9c")

  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2_empty_col.csv"),
                      header = T)
  exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
  expect_known_hash(exp_type_2, "f8d57b6e9c")

  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2_empty_row.csv"),
                      header = T)
  exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
  expect_known_hash(exp_type_2, "516b7aee57")

  # calibration data
  cal_type_2 <- fread(paste0(prefix, "testdata/cal_type_2.csv"))
  cal_type_2 <- clean_dt(cal_type_2, "calibration", 2, logfilename)[["dat"]]
  expect_known_hash(cal_type_2, "d32a53505b")
})

test_that(
  desc = "wrong description",
  code = {
  # type 1 data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
  expect_error(clean_dt(cal_type_1, "calibraRAtion", 1, logfilename))

  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  expect_error(clean_dt(exp_type_1, "experiRINKLmental", 1, logfilename))

  # type 2 data
  cal_type_2 <- fread(paste0(prefix, "testdata/cal_type_2.csv"))
  expect_error(clean_dt(cal_type_2, "calibraRAtion", 2, logfilename))

  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2.csv"))
  expect_error(clean_dt(exp_type_2, "experiRINKLmental", 2, logfilename))
})

# wrong type
test_that(
  desc = "wrong type",
  code = {
  # type 1 data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1.csv"))
  expect_null(clean_dt(cal_type_1, "calibration", 3, logfilename))
  expect_null(clean_dt(cal_type_1, "calibration", "a", logfilename))

  exp_type_1 <- fread(paste0(prefix, "testdata/exp_type_1.csv"))
  expect_null(clean_dt(exp_type_1, "experimental", 65, logfilename))
  expect_null(clean_dt(exp_type_1, "experimental", "tre", logfilename))

  # type 2 data
  cal_type_2 <- fread(paste0(prefix, "testdata/cal_type_2.csv"))
  expect_null(clean_dt(cal_type_2, "calibration", 3, logfilename))
  expect_null(clean_dt(cal_type_2, "calibration", "a", logfilename))

  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2.csv"))
  expect_null(clean_dt(exp_type_2, "experimental", 65, logfilename))
  expect_null(clean_dt(exp_type_2, "experimental", "tre", logfilename))
})

# wrong first col
test_that(
  desc = "wrong first column in calibration data type 1",
  code = {
  # type 1 data
  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1_wrong_col_1.csv"))
  expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))

  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1_wrong_col_1_2.csv"))
  expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))

  cal_type_1 <- fread(paste0(prefix, "testdata/cal_type_1_less4.csv"))
  expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))
})

# heterogenous cpg-sites per locus
test_that(
  desc = "heterogenous cpg-sites per locus in type 2 data",
  code = {
  cal_type_2 <- fread(paste0(prefix, "testdata/cal_type_2_heterogenous.csv"))
  expect_null(clean_dt(cal_type_2, "calibration", 2, logfilename))

  exp_type_2 <- fread(paste0(prefix, "testdata/exp_type_2_heterogenous.csv"))
  expect_null(clean_dt(exp_type_2, "experimental", 2, logfilename))

  expect_true(file.remove(paste0(prefix, "log.txt")))
})

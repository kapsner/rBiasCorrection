context("dataimport filechecks")

prefix <- tempdir()

# the writeLog-function needs the logfilename
logfilename <- paste0(prefix, "/log.txt")

library(data.table)

test_that(
  desc = "test normal function of file import of type 1",
  code = {

    local_edition(3)
    # experimental data
    exp_type_1 <- fread("./testdata/exp_type_1.csv")
    exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)
    expect_type(exp_type_1, "list")
    expect_snapshot_value(
      x = table_prep(exp_type_1[["dat"]]),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_1[["dat"]], "ab5287b084")
    #   expect_known_hash(exp_type_1[["dat"]], "8c93f2d28f")
    # }, class = "error", regexp = "ab5287b084|8c93f2d28f")
    # nolint end

    exp_type_1 <- fread("./testdata/exp_type_1_empty_col.csv", header = T)
    exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(exp_type_1),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_1, "ab5287b084")
    #   expect_known_hash(exp_type_1, "8c93f2d28f")
    # }, class = "error", regexp = "ab5287b084|8c93f2d28f")
    # nolint end

    exp_type_1 <- fread("./testdata/exp_type_1_empty_row.csv")
    exp_type_1 <- clean_dt(exp_type_1, "experimental", 1, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(exp_type_1),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_1, "29fa13cb24")
    #   expect_known_hash(exp_type_1, "1fb62d8498")
    # }, class = "error", regexp = "29fa13cb24|1fb62d8498")
    # nolint end

    # calibration data
    cal_type_1 <- fread("./testdata/cal_type_1.csv")
    cal_type_1 <- clean_dt(cal_type_1, "calibration", 1, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(cal_type_1),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(cal_type_1, "869cda3040")
    #   expect_known_hash(cal_type_1, "23f21fc354")
    # }, class = "error", regexp = "869cda3040|23f21fc354")
    # nolint end
  })

test_that(
  desc = "test normal function of file import of type 2",
  code = {

    local_edition(3)
    # experimental data
    exp_type_2 <- fread("./testdata/exp_type_2.csv")
    exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(exp_type_2),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_2, "f8d57b6e9c")
    #   expect_known_hash(exp_type_2, "78b933c115")
    # }, class = "error", regexp = "f8d57b6e9c|78b933c115")
    # nolint end

    exp_type_2 <- fread("./testdata/exp_type_2_empty_col.csv", header = T)
    exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(exp_type_2),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_2, "f8d57b6e9c")
    #   expect_known_hash(exp_type_2, "78b933c115")
    # }, class = "error", regexp = "f8d57b6e9c|78b933c115")
    # nolint end

    exp_type_2 <- fread("./testdata/exp_type_2_empty_row.csv")
    exp_type_2 <- clean_dt(exp_type_2, "experimental", 2, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(exp_type_2),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(exp_type_2, "eba558a583")
    #   expect_known_hash(exp_type_2, "516b7aee57")
    # }, class = "error", regexp = "eba558a583|516b7aee57")
    # nolint end

    # calibration data
    cal_type_2 <- fread("./testdata/cal_type_2.csv")
    cal_type_2 <- clean_dt(cal_type_2, "calibration", 2, logfilename)[["dat"]]
    expect_snapshot_value(
      x = table_prep(cal_type_2),
      style = "json2",
      cran = FALSE,
      tolerance = 10e-3,
      ignore_function_env = TRUE
    )
    # nolint start
    # expect_error({
    #   expect_known_hash(cal_type_2, "296406ae1c")
    #   expect_known_hash(cal_type_2, "d32a53505b")
    # }, class = "error", regexp = "296406ae1c|d32a53505b")
    # nolint end
  })

test_that(
  desc = "wrong description",
  code = {
    # type 1 data
    cal_type_1 <- fread("./testdata/cal_type_1.csv")
    expect_error(clean_dt(cal_type_1, "calibraRAtion", 1, logfilename))

    exp_type_1 <- fread("./testdata/exp_type_1.csv")
    expect_error(clean_dt(exp_type_1, "experiRINKLmental", 1, logfilename))

    # type 2 data
    cal_type_2 <- fread("./testdata/cal_type_2.csv")
    expect_error(clean_dt(cal_type_2, "calibraRAtion", 2, logfilename))

    exp_type_2 <- fread("./testdata/exp_type_2.csv")
    expect_error(clean_dt(exp_type_2, "experiRINKLmental", 2, logfilename))
  })

# wrong type
test_that(
  desc = "wrong type",
  code = {
    # type 1 data
    cal_type_1 <- fread("./testdata/cal_type_1.csv")
    expect_null(clean_dt(cal_type_1, "calibration", 3, logfilename))
    expect_null(clean_dt(cal_type_1, "calibration", "a", logfilename))

    exp_type_1 <- fread("./testdata/exp_type_1.csv")
    expect_null(clean_dt(exp_type_1, "experimental", 65, logfilename))
    expect_null(clean_dt(exp_type_1, "experimental", "tre", logfilename))

    # type 2 data
    cal_type_2 <- fread("./testdata/cal_type_2.csv")
    expect_null(clean_dt(cal_type_2, "calibration", 3, logfilename))
    expect_null(clean_dt(cal_type_2, "calibration", "a", logfilename))

    exp_type_2 <- fread("./testdata/exp_type_2.csv")
    expect_null(clean_dt(exp_type_2, "experimental", 65, logfilename))
    expect_null(clean_dt(exp_type_2, "experimental", "tre", logfilename))
  })

# wrong first col
test_that(
  desc = "wrong first column in calibration data type 1",
  code = {
    # type 1 data
    cal_type_1 <- fread("./testdata/cal_type_1_wrong_col_1.csv")
    expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))

    cal_type_1 <- fread("./testdata/cal_type_1_wrong_col_1_2.csv")
    expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))

    cal_type_1 <- fread("./testdata/cal_type_1_less4.csv")
    expect_null(clean_dt(cal_type_1, "calibration", 1, logfilename))
  })

# heterogenous cpg-sites per locus
test_that(
  desc = "heterogenous cpg-sites per locus in type 2 data",
  code = {
    cal_type_2 <- fread("./testdata/cal_type_2_heterogenous.csv")
    expect_null(clean_dt(cal_type_2, "calibration", 2, logfilename))

    exp_type_2 <- fread("./testdata/exp_type_2_heterogenous.csv")
    expect_null(clean_dt(exp_type_2, "experimental", 2, logfilename))

    expect_true(file.remove(paste0(prefix, "/log.txt")))
  })

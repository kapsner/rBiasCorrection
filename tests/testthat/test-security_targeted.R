library(testthat)

test_that("Path validation (.is_path_safe) works correctly", {
  # We need to access the private function
  is_path_safe <- rBiasCorrection:::.is_path_safe

  # Safe paths
  expect_true(is_path_safe(tempdir()))
  expect_true(is_path_safe(getwd()))
  expect_true(is_path_safe(file.path(tempdir(), "subdir")))
  expect_true(is_path_safe(file.path(getwd(), "subdir")))

  # Unsafe paths
  expect_false(is_path_safe("/"))
  if (.Platform$OS.type == "windows") {
    expect_false(is_path_safe("C:\\"))
  }
  expect_false(is_path_safe("/etc"))
  expect_false(is_path_safe("/home"))

  # Path traversal
  expect_false(is_path_safe(file.path(getwd(), "../../../etc/passwd")))

  # Empty/NULL
  expect_false(is_path_safe(""))
  expect_false(is_path_safe(NULL))
})

test_that("biascorrection stops on unsafe paths", {
  # Mocking data files
  exp_file <- tempfile(fileext = ".csv")
  cal_file <- tempfile(fileext = ".csv")
  data.table::fwrite(rBiasCorrection::example.data_experimental$dat, exp_file)
  data.table::fwrite(rBiasCorrection::example.data_calibration$dat, cal_file)

  # Unsafe plotdir
  expect_error(
    biascorrection(
      experimental = exp_file,
      calibration = cal_file,
      samplelocusname = "Test",
      plotdir = "/",
      parallel = FALSE
    ),
    "Security risk: `plotdir=/` is outside of allowed paths"
  )

  # Unsafe csvdir
  expect_error(
    biascorrection(
      experimental = exp_file,
      calibration = cal_file,
      samplelocusname = "Test",
      csvdir = "/etc",
      parallel = FALSE
    ),
    "Security risk: `csvdir=/etc/` is outside of allowed paths"
  )

  # Unsafe logfilename
  expect_error(
    biascorrection(
      experimental = exp_file,
      calibration = cal_file,
      samplelocusname = "Test",
      logfilename = "/var/log/system.log",
      parallel = FALSE
    ),
    "Security risk: `logfilename=/var/log/system.log` is outside of allowed paths"
  )
})

test_expect_clean_up_safety <- function() {
  # This is hard to test directly without actually trying to delete something sensitive,
  # which we shouldn't do. But we can test that it doesn't fail on safe paths.
  plotdir <- file.path(tempdir(), "test_plotdir")
  csvdir <- file.path(tempdir(), "test_csvdir")
  dir.create(plotdir)
  dir.create(csvdir)

  expect_silent(rBiasCorrection::clean_up(plotdir, csvdir))
  expect_false(dir.exists(plotdir))
  expect_false(dir.exists(csvdir))
}

test_that("clean_up handles safe paths correctly", {
  test_expect_clean_up_safety()
})

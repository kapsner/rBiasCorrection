context("lints")

prefix <- "./"

#% prefix <- "../../"

test_that(
  desc = "test lints",
  code = {
    lintlist <- list(
      "R" = list(
        "better_model.R" = NULL,
        "biascorrection.R" = NULL,
        "clean_dt.R" = "cyclomatic complexity",
        "cubic.R" = NULL,
        "hyperbolic.R" = NULL,
        "plotting.R" = NULL,
        "regression.R" = NULL,
        "regression_solving.R" = "cyclomatic complexity",
        "statistics_list.R" = NULL,
        "type2files.R" = NULL,
        "utils.R" = NULL
      ),
      "tests/testthat" = list(
        "test-algorithm_minmax_FALSE.R" = NULL,
        "test-algorithm_minmax_TRUE.R" = NULL,
        "test-biascorrection.R" = NULL,
        "test-clean_dt.R" = NULL,
        "test-create_aggregated.R" = NULL,
        "test-hyperbolic.R" = NULL,
        "test-lints.R" = NULL,
        "test-utils.R" = NULL
      )
    )
    for (directory in names(lintlist)) {
      print(directory)
      for (fname in names(lintlist[[directory]])) {
        print(fname)
        lintr::expect_lint(
          file = paste0(prefix,
                        directory,
                        "/",
                        fname),
          checks = lintlist[[directory]][[fname]])
      }
    }
  })

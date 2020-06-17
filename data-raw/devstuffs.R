packagename <- "rBiasCorrection"

# remove existing description object
unlink("DESCRIPTION")
# Create a new description object
my_desc <- desc::description$new("!new")
# Set your package name
my_desc$set("Package", packagename)
#Set your name
my_desc$set_authors(c(
  person("Lorenz A.", "Kapsner", email = "lorenz.kapsner@gmail.com", role = c('cre', 'aut'),
         comment = c(ORCID = "0000-0003-1866-860X")),
  person("Evgeny A.", "Moskalev", role = "aut")))
# Remove some author fields
my_desc$del("Maintainer")
# Set the version
my_desc$set_version("0.1.7")
# The title of your package
my_desc$set(Title = "A Package to Correct Bias in DNA Methylation Analyses")
# The description of your package
my_desc$set(Description = "Correct Bias in Quantitative DNA Methylation Analyses.")
# The description of your package
my_desc$set("Date/Publication" = paste(as.character(Sys.time()), "UTC"))
# The urls
my_desc$set("URL", "https://github.com/kapsner/rBiasCorrection")
my_desc$set("BugReports",
            "https://github.com/kapsner/rBiasCorrection/issues")

# Vignette Builder
my_desc$set("VignetteBuilder" = "knitr")

# License
my_desc$set("License", "GPL-3")
# Save everyting
my_desc$write(file = "DESCRIPTION")

# License
#usethis::use_gpl3_license(name="Lorenz Kapsner")

# Depends
usethis::use_package("R", min_version = "2.10", type="Depends")

# Imports
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html
usethis::use_package("data.table", type="Imports")
usethis::use_package("shiny", type="Imports")
usethis::use_package("ggplot2", type="Imports")
usethis::use_package("magrittr", type="Imports")
usethis::use_package("polynom", type="Imports")
usethis::use_package("nls2", type="Imports")
usethis::use_package("ggpubr", type="Imports")
usethis::use_package("stats", type="Imports")

# Suggests
usethis::use_package("testthat", type = "Suggests")
usethis::use_package("lintr", type = "Suggests")
usethis::use_package("knitr", type = "Suggests")


# gitignore
usethis::use_git_ignore("/*")
usethis::use_git_ignore("/*/")
usethis::use_git_ignore("*.log")
usethis::use_git_ignore("!/.gitignore")
usethis::use_git_ignore("!/.Rbuildignore")
usethis::use_git_ignore("!/.gitlab-ci.yml")
usethis::use_git_ignore("!/data-raw/")
usethis::use_git_ignore("!/data/")
usethis::use_git_ignore("!/inst/")
usethis::use_git_ignore("/data-raw/*")
usethis::use_git_ignore("!/data-raw/devstuffs.R")
usethis::use_git_ignore("!/DESCRIPTION")
usethis::use_git_ignore("!/FAQ.md")
usethis::use_git_ignore("!/LICENSE.md")
usethis::use_git_ignore("!/man/")
usethis::use_git_ignore("!NAMESPACE")
usethis::use_git_ignore("!/R/")
usethis::use_git_ignore("!/README.md")
usethis::use_git_ignore("!/vignettes/")
usethis::use_git_ignore("!/tests/")
usethis::use_git_ignore("/tests/*")
usethis::use_git_ignore("!/tests/testthat/")
usethis::use_git_ignore("/tests/testthat/*")
usethis::use_git_ignore("!/tests/testthat/test-*.R")
usethis::use_git_ignore("!/tests/testthat.R")
usethis::use_git_ignore("/tests/testthat/testdata/")
usethis::use_git_ignore("/.Rhistory")
usethis::use_git_ignore("!/*.Rproj")
usethis::use_git_ignore("/.Rproj*")
usethis::use_git_ignore("/.RData")
usethis::use_git_ignore("/.vscode")
usethis::use_git_ignore("!/.lintr")
usethis::use_git_ignore("!/.github/")
usethis::use_git_ignore("!/tic.R")

usethis::use_build_ignore("csvdir")
usethis::use_build_ignore("plotdir")
usethis::use_build_ignore("log.txt")
usethis::use_build_ignore("revdep")
usethis::use_build_ignore("FAQ.md")
usethis::use_build_ignore("FAQ.html")
usethis::use_build_ignore("tests/testthat/csvdir")
usethis::use_build_ignore("tests/testthat/plotdir")
usethis::use_build_ignore(".vscode")
usethis::use_build_ignore(".lintr")
usethis::use_build_ignore("tic.R")
usethis::use_build_ignore(".github")


# experimental = "../19_PCR-bias/data/example_data/type1/example_data_type1_experimentaldata.csv"
# calibration = "../19_PCR-bias/data/example_data/type1/example_data_type1_calibrationdata.csv"
# data
example.data_experimental <- data.table::fread("./tests/testthat/testdata/exp_type_1.csv")
example.data_calibration <- data.table::fread("./tests/testthat/testdata/cal_type_1.csv")

# vec <- colnames(example.data_experimental)
# vec <- vec[2:length(vec)]
# example.data_experimental[, (vec) := lapply(.SD, function(x) {
#   as.numeric(gsub("\\,", ".", x))
#   }), .SDcols=vec]
#
# data.table::fwrite(example.data_experimental, "./tests/testthat/testdata/exp_type_1.csv")
#
# vec <- colnames(example.data_calibration)
# example.data_calibration[, (vec) := lapply(.SD, function(x) {
#   as.numeric(gsub("\\,", ".", x))
# }), .SDcols=vec]
#
# data.table::fwrite(example.data_calibration, "./tests/testthat/testdata/cal_type_1.csv")

example.data_experimental <- rBiasCorrection::clean_dt(example.data_experimental,
                                                       "experimental",
                                                       1,
                                                       "./logfile.txt")
example.data_calibration <- rBiasCorrection::clean_dt(example.data_calibration,
                                                      "calibration",
                                                      1,
                                                      "./logfile.txt")
file.remove("./logfile.txt")
usethis::use_data(example.data_experimental, example.data_calibration, internal = F)


# BiasCorrection(experimental = "../19_PCR-bias/data/example_data/type1/example_data_type1_experimentaldata.csv", calibration = "../19_PCR-bias/data/example_data/type1/example_data_type1_calibrationdata.csv", samplelocusname = "Test")
#covr::package_coverage()
#lintr::lint_package()
#styler::style_pkg()


# TODO solve minmax cubic eq

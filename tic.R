# installs dependencies, runs R CMD check, runs covr::codecov()
do_package_checks()

get_stage("script") %>%
    add_code_step(install.packages("rmarkdown"))

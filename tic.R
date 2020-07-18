# installs dependencies, runs R CMD check, runs covr::codecov()
do_package_checks()

get_stage("install") %>%
    add_code_step(remotes::update_packages(packages = "rmarkdown", upgrade = "always"))

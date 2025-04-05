# https://github.com/Rdatatable/data.table/issues/5658
Sys.setenv("OMP_THREAD_LIMIT" = 2)
Sys.setenv("Ncpu" = 2)

library(testthat)
library(rBiasCorrection)

local_edition(3)

test_check("rBiasCorrection")

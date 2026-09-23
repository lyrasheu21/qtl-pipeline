# Run the unit tests for the shared R helpers:  Rscript tests/run_tests.R
library(testthat)
suppressPackageStartupMessages(library(data.table))
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)

# Run the unit tests without installing the package:
#   Rscript qtlpost/tests/run_tests.R      (from the repository root)
library(testthat)
here <- if (dir.exists("qtlpost/tests/testthat")) "qtlpost/tests/testthat" else "tests/testthat"
test_dir(here, reporter = "summary", stop_on_failure = TRUE)

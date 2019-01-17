# Test
library(testthat) 

# includes go here
#source("*.R")

#test_results <- test_dir("tests/testthat", reporter="summary")
test_results <- test_dir("tests/testthat", reporter = LocationReporter)
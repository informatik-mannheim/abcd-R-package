# Test
# Not sure, if this is allowed/intended:
#source("../../R/skewCalc.R")

# Single skew calculations go here:

test_that("Test Partition A", {
  res = subseqpart("A", 1, 1)
  expect_equal(res, "A")
})

test_that("Test Partition A default values", {
  res = subseqpart("A")
  expect_equal(res, "A")
})

test_that("Test Partition AT", {
  seq = "AT"
  res = subseqpart(seq, 2, 1)
  expect_equal(res, "A")
  res = subseqpart(seq, 2, 2)
  expect_equal(res, "T")
})

test_that("Test Partition ATC", {
  seq = "ATC"
  res = subseqpart(seq, 2, 1)
  expect_equal(res, "A")
  res = subseqpart(seq, 2, 2)
  expect_equal(res, "T")
})

test_that("Test Partition ATCG", {
  seq = "ATCG"
  res = subseqpart(seq, 2, 1)
  expect_equal(res, "AT")
  res = subseqpart(seq, 2, 2)
  expect_equal(res, "CG")
})

test_that("Test Partition ATCGA", {
  seq = "ATCGA"
  res = subseqpart(seq, 2, 1)
  expect_equal(res, "AT")
  res = subseqpart(seq, 2, 2)
  expect_equal(res, "CG")
})
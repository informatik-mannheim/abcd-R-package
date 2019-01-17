# Test
# Not sure, if this is allowed/intended:
# source("../../R/skewCalc.R")

# Single skew calculations go here:

test_that("Test Skew at", {
  seq <- "at"
  res <- skew(seq)
  expect_equal(res, c(at=0, cg=NaN))
})

test_that("Test Skew cg", {
  seq <- "cg"
  res <- skew(seq)
  expect_equal(res, c(at=NaN, cg=0))
})

test_that("Test Skew atcg", {
  seq <- "atcg"
  res <- skew(seq)
  expect_equal(res, c(at=0, cg=0))
})

test_that("Test Skew atcga", {
  seq <- "atcga"
  res <- skew(seq)
  expect_equal(res, c(at=1/3, cg=0), tolerance = 0.001)
})

# n-plets
test_that("Test 2-plet Skew atcatc", {
  # k = 1: a c t   
  # -> AT = (1-1)/(1+1) = 0; CG = (1-0)/(1+0) = 1
  # k = 2:  t a c
  # -> AT = (1-1)/(1+1) = 0; CG = (1-0)/(1+0) = 1
  seq <- "atcatc"
  res <- nplet.skew(seq, 2)
  skewsAT <- res$at
  skewsCG <- res$cg
  expect_equal(length(skewsAT), 2)
  expect_equal(length(skewsCG), 2)
  expect_equal(skewsAT, c("1"=0, "2"=0))
  expect_equal(skewsCG, c("1"=1, "2"=1))
})

test_that("Test 3-plet Skew atgatccgcgatgta", {
  # k = 1: a  a  c  c  g
  # -> AT = (2-0)/(2+0) = 1; CG = (2-1)/(2+1) = 1/3
  # k = 2:  t  t  g  a  t
  # -> AT = (1-3)/(1+3) = -1/2; CG = (0-1)/0+1) = -1
  # k = 3:   g  c  c  t  a
  # -> AT = (1-1)/(1+1) = 0; CG = (2-1)/2+1) = 1/3
  seq <- "atgatccgccatgta"
  res <- nplet.skew(seq, 3)
  skewsAT <- res$at
  skewsCG <- res$cg
  expect_equal(length(skewsAT), 3)
  expect_equal(length(skewsCG), 3)
  expect_equal(skewsAT, c("1"=1, "2"=-1/2, "3"=0))
  expect_equal(skewsCG, c("1"=1/3, "2"=-1, "3"=1/3))
})
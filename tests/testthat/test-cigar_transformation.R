library(testthat)
library(mappedreadalignment)

example_cigar <- "3H3M2D2M"
result <- cigar_transform(example_cigar)

test_that("cigar_transform correct CIGAR transformation", {
  expect_equal(result,c("H","H","H","M","M","M","D","D","M","M") )})

example_wrong_cigar <- "3H3M2J2K"

test_that("cigar_transform raise an error for invalid CIGAR operations", {
  expect_error(cigar_transform(example_wrong_cigar))
  "Invalid operations in CIGAR string: J, K"})

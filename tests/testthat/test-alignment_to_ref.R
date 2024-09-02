library(testthat)
library(mappedreadalignment)

example_read <- list(rname = "chr6",
                     pos = 32401000, cigar = "2S5M2D4M",
                     seq = "AAAGATCGACC", qname = "read4")

result <- alignment_to_ref(example_read)


test_that("alignment_to_ref correct Read_ID", {
  expect_equal(result$Read_ID, "read4")})

test_that("alignment_to_ref correct CIGAR", {
  expect_equal(result$CIGAR, "2S5M2D4M")})

test_that("alignment_to_ref correct Reference_genome", {
  expect_equal(result$Alignment$Reference_genome, "  AGATCGAGACC")})

test_that("alignment_to_ref correct Query_read", {
  expect_equal(result$Alignment$Query_read, "AAAGATC--GACC")})



example_wrong_read_cigar <- list(rname = "chr6",
                                 pos = 32401000, cigar = "2S5M2K4J",
                                 seq = "AAAGATCGACC", qname = "read4")

test_that("alignment_to_ref raise an error for invalid CIGAR operations", {
  expect_error(alignment_to_ref(example_wrong_read_cigar))
  "Invalid operations in CIGAR string: K, J"})

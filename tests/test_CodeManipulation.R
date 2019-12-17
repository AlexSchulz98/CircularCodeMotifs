source("CodeManipulation.R")

library(testthat)
library(Biostrings)
library(ccmotif)

test_that(
  "part of circular function with true values", 
  {
  ranCode = sample(1:216,1)
  setX = codes.c3[[ranCode]]
  ispartof = unlist(setX[2])

  ranCodon = sample(1:length(ispartof),1)

  expect_true(partOfCircularCode(ispartof[[ranCodon]],codes.c3[[ranCode]]))
})

test_that(
  "part of circular function with false value", 
  {
    ranCode = sample(1:216,1)
    
    expect_false(partOfCircularCode("AAA",codes.c3[[ranCode]]))
  })
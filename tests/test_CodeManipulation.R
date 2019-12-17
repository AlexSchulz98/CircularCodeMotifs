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

test_that(
  "codon subsitution function - 1 max value - take max value", 
  {
    #AAG = 9, AAT = 8, AAC = 8
    expect_equal(codonSubstitution(oldCodon = "AAA", newCodons = c("AAG","AAT","AAC")),"AAG")
  })

test_that(
  "codon subsitution function - 2 max values - take first max value", 
  {
    #AAT = 8, AAC = 8, AGA = 7
    expect_equal(codonSubstitution(oldCodon = "AAA", newCodons = c("AGA","AAC","AAT")),"AAC")
  })

test_that(
  "codon subsitution function - zero matches - take first entry", 
  {
    #CTC = 0, CCT = 0, TCT = 0
    expect_equal(codonSubstitution(oldCodon = "AAA", newCodons = c("CTC","CCT","TCT")),"CTC")
  })

test_that(
  "amine change function - A", 
  {
    expect_equal(amineChange("A","G"),1)
    expect_equal(amineChange("A","T"),0)
    expect_equal(amineChange("A","C"),0)
  })

test_that(
  "amine change function - G", 
  {
    expect_equal(amineChange("G","A"),1)
    expect_equal(amineChange("G","T"),0)
    expect_equal(amineChange("G","C"),0)
  })

test_that(
  "amine change function - C", 
  {
    expect_equal(amineChange("C","G"),0)
    expect_equal(amineChange("C","A"),0)
    expect_equal(amineChange("C","T"),1)
  })

test_that(
  "amine change function - T", 
  {
    expect_equal(amineChange("T","G"),0)
    expect_equal(amineChange("T","A"),0)
    expect_equal(amineChange("T","C"),1)
  })

test_that(
  "amino acid substitution function", 
  {
    expect_equal(amineChange("T","G"),0)
    expect_equal(amineChange("T","A"),0)
    expect_equal(amineChange("T","C"),1)
  })





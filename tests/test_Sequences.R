source("/Projekte/BA Circular Code/Scripts/Sequences.R")

library(testthat)
library(Biostrings)
library(ccmotif)


test_that(
  "change reading frame function - neagtive value",
  {
    seq = DNAString("AAATTTCCCGGG")
    
    expect_error(changeReadingFrame(-1, seq))
  })

test_that(
  "change reading frame function - by 0",
  {
    seq = DNAString("AAATTTCCCGGG")
    
    expected = "AAATTTCCCGGG"
    actual = toString(changeReadingFrame(0, seq))
    
    expect_equal(expected, actual)
  })

test_that(
  "change reading frame function - by 1",
  {
    seq = DNAString("AAATTTCCCGGG")
    
    expected = "AATTTCCCGGG"
    actual = toString(changeReadingFrame(1, seq))
    
    expect_equal(expected, actual)
  })

test_that(
  "change reading frame function - by 2",
  {
    seq = DNAString("AAATTTCCCGGG")
    
    expected = "ATTTCCCGGG"
    actual = toString(changeReadingFrame(2, seq))
    
    expect_equal(expected, actual)
  })

test_that(
  "change reading frame function - by 3",
  {
    seq = DNAString("AAATTTCCCGGG")
    
    expected = "TTTCCCGGG"
    actual = toString(changeReadingFrame(3, seq))
    
    expect_equal(expected, actual)
  })


test_that(
  "get IUPAC Codes Count function - no IUPAC in sequence",
  {
    seqSet = DNAStringSet("AAATTTCCCGGG")
    
    expected = 0
    actual = getIUPACCodesCount(seqSet, 1, 1)
    
    expect_equal(expected, actual)
  })

test_that(
  "get IUPAC Codes Count function - IUPAC in sequence",
  {
    seqA = DNAString("AAATTTGGGCCC")
    seqB = DNAString("KKKTTTGGGYYY")
    seqSet = DNAStringSet(x=c(seqA,seqB),start = c(1,13), width = 12)
    
    expected = 1
    actual = getIUPACCodesCount(seqSet, 1, 2)
    
    expect_equal(expected, actual)
  })


test_that(
  "delete IUPAC Sequences function - no deletions",
  {
    seqA = DNAString("AAATTTGGGCCC")
    seqB = DNAString("CCCGGGTTTAAA")
    seqSet = DNAStringSet(x=c(seqA,seqB),start = c(1,13), width = 12)
    
    expected = 2
    actual = deleteIUPACSequences(seqSet)
    actual = length(actual)
    
    expect_equal(expected, actual)
  })

test_that(
  "delete IUPAC Sequences function - delete 1 sequence",
  {
    seqA = DNAString("AAATTTGGGCCC")
    seqB = DNAString("KKKTTTGGGYYY")
    seqSet = DNAStringSet(x=c(seqA,seqB),start = c(1,13), width = 12)
    
    expected = 1
    actual = deleteIUPACSequences(seqSet)
    actual = length(actual)
    
    expect_equal(expected, actual)
  })


test_that(
  "shift code function - shift by 1",
  {
    code = 23
    expected = c("ACA", "ATA", "CCA", "TCA" ,"TTA", "AGC" ,"TCC", "TGC" ,"AAG", "ACG", "AGG", "ATG", "CCG" ,"GCG", "GTG", "TAG", "TCG", "TTG", "ACT" ,"TCT")
    
    actual = shiftCode(code, 1)
    
    expect_equal(expected, actual)
  })

test_that(
  "shift code function - shift by 2",
  {
    code = 23
    expected = c("CAA", "TAA" ,"CAC" ,"CAT", "TAT" ,"GCA", "CCT", "GCT", "AGA" ,"CGA", "GGA","TGA" ,"CGC" ,"CGG", "TGG", "AGT" ,"CGT" ,"TGT" ,"CTA", "CTT")
    
    actual = shiftCode(code, 2)
    
    expect_equal(expected, actual)
  })


test_that(
  "part of circular code list - no codons part of the circular code",
  {
    codons = c("AAA","TTT","GGG")
    code = 23
    
    expected = c(F,F,F)
    actual = partOfCircularCode_List(codons, code)
    
    expect_equal(expected, actual)
  })

test_that(
  "part of circular code list - with 2 / 4 codons of circular code",
  {
    codons = c("AAA","TAC","GGG", "TTC")
    code = 23
    
    expected = c(F,T,F,T)
    actual = partOfCircularCode_List(codons, code)
    
    expect_equal(expected, actual)
  })
















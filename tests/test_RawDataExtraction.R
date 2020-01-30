source("/Projekte/BA Circular Code/RawDataExtraction.R")
source("/Projekte/BA Circular Code/CodeManipulation.R")

library(testthat)
library(Biostrings)
library(ccmotif)

# 64 codons
CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
# 20 amino acids + 1 *  for stop codons
AMINOACIDS = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

test_that(
  "generate empty table function",
  {
    dim = sample(1:100,1)
    
    ar = generateEmptyTable(dim, dim, 1:dim)
    
    expected = dim*dim
    actual = length(ar)
    
    expect_equal(expected, actual)
    expect_equal(sum(ar), 0)
  })

test_that(
  "codon count function",
  {
    ar = generateEmptyTable(64, 64, CODONS)
    seqA = DNAString("AAATTT")
    seqB = DNAString("TTTTTT")
    
    expected = ar
    expected["AAA","TTT"] = 1
    expected["TTT","TTT"] = 1
    
    actual = codonCount(ar,seqA,seqB)
    
    expect_equal(expected, actual)
  })

test_that(
  "amino count function",
  {
    ar = generateEmptyTable(21, 21, AMINOACIDS)
    seqA = DNAString("AAATTT")
    seqB = DNAString("TTTTTT")
    
    expected = ar
    expected["K","F"] = 1
    expected["F","F"] = 1
    
    actual = aminoCount(ar,seqA,seqB)
    
    expect_equal(expected, actual)
  })

test_that(
  "unchanged non CC codons function",
  {
    code = codes.c3[[23]]
    codes = code[[2]]
    ar = generateEmptyTable(64, 64, CODONS)
    ar[toString(codes[1]),toString(codes[1])] = 30
    ar["TTT","TTT"] = 20
    ar["CCC","AAA"] = 10
    
    expected = 20
    
    actual = unchangednonCCCodons(ar, code)
    
    expect_equal(expected, actual)
  })


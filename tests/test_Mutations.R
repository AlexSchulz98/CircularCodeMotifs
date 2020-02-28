source("/Projekte/BA Circular Code/Scripts/Mutations.R")

library(testthat)
library(Biostrings)
library(ccmotif)

test_that(
  "part of circular function - true",
  {
  ranCode = sample(1:216,1)
  setX = codes.c3[[ranCode]]
  ispartof = unlist(setX[2])
  ranCodon = sample(1:length(ispartof),1)
  
  actual = partOfCircularCode(ispartof[[ranCodon]],codes.c3[[ranCode]])

  expect_true(actual)
})

test_that(
  "part of circular function - false",
  {
    ranCode = sample(1:216,1)
    
    actual = partOfCircularCode("AAA",codes.c3[[ranCode]])

    expect_false(actual)
  })


test_that(
  "get codes for AA function - normal codon",
  {
    aa = AAString("L")
    
    expected = c("TTA","TTG","CTA","CTG","CTT","CTC")
    actual = getCodesForAA(aa)
    
    expect_equal(expected, actual)
  })

test_that(
  "get codes for AA function - stop codon",
  {
    aa = AAString("*")
    
    expected = c("TAA","TAG","TGA")
    actual = getCodesForAA(aa)
    
    expect_equal(expected, actual)
  })


test_that(
  "get circular codes function",
  {
    setX = codes.c3[[23]]
    codons = c("ATT","GGG","AAA","CTC")
    
    expected = c("ATT","CTC")
    actual = getCircularCodes(codons, setX)
    
    expect_equal(expected, actual)
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
  "codon subsitution function - 1 max value - take max value",
  {
    #Scores: AAG = 9, AAT = 8, AAC = 8
    expected = "AAG"
    actual = codonSubstitution(oldCodon = "AAA", newCodons = c("AAG","AAT","AAC"))
    
    expect_equal(expected, actual)
  })

test_that(
  "codon subsitution function - 2 max values - take first max value",
  {
    #Scores: AAT = 8, AAC = 8, AGA = 7
    expected = "AAC"
    actual = codonSubstitution(oldCodon = "AAA", newCodons = c("AGA","AAC","AAT"))
    
    expect_equal(expected, actual)
  })

test_that(
  "codon subsitution function - zero matches - take first entry",
  {
    #Scores: CTC = 0, CCT = 0, TCT = 0
    expected = "CTC"
    actual = codonSubstitution(oldCodon = "AAA", newCodons = c("CTC","CCT","TCT"))
    
    expect_equal(expected, actual)
  })


test_that(
  "get amino acids coded by X function",
  {
    setX = codes.c3[[23]]
    
    expected = toString(AAStringSet(c("N","T","I","Q","L","M","E","D","A","G","V","Y","F")))
    actual = toString(getAminoAcidsCodedByX(setX))
    
    expect_equal(expected, actual)
  })


test_that(
  "amino acid substitution function - low threshold for change",
  {
    expected = "E"
    actual = toString(aminoAcidSubstitution(codon = DNAString("AAA"),setX = codes.c3[[26]],threshold = 0))
    
    expect_equal(expected, actual)
  })

test_that(
  "amino acid substitution function - high treshold for no change",
  {
    expected = "K"
    actual = toString(aminoAcidSubstitution(codon = DNAString("AAA"),setX = codes.c3[[26]],threshold = 2))
    
    expect_equal(expected, actual)
  })


test_that(
  "mutate sequence function - without amino acid substitution",
  {
    seq = DNAString("AAATTTGGGCCC")
    expected = DNAString("AAATTCGGTCCC")
    actual = mutateSequence(seq, codes.c3[[23]])
    
    expect_equal(expected, actual)
  })

test_that(
  "mutate sequence function - with amino acid substitution - too high threshold",
  {
    seq = DNAString("CGG")
    expected = DNAString("CGG")
    actual = mutateSequence(seq, codes.c3[[23]], subAminos = TRUE, threshold = 5)
    
    expect_equal(expected, actual)
  })

test_that(
  "mutate sequence function - with amino acid substitution - threshold zero",
  {
    seq = DNAString("CGG")
    expected = DNAString("CGA")
    actual = mutateSequence(seq, codes.c3[[23]], subAminos = TRUE, threshold = 0)
    
    expect_equal(expected, actual)
  })

test_that(
  "mutate sequence function - with amino acid substitution - low threshold, still take highest score",
  {
    seq = DNAString("CGG")
    expected = DNAString("CGA")
    actual = mutateSequence(seq, codes.c3[[23]], subAminos = TRUE, threshold = -4)
    
    expect_equal(expected, actual)
  })

























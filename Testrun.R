
library(Biostrings)
library(ccmotif)
source("ChangeCodons.R")

seqSet = readDNAStringSet("cds/ena-sars.fasta") # RNA sequence set
seq=seqSet[[11]] #Single RNA sequence
codons = codons(seq) # List of codons in rnaseq

setX = codes.c3[[26]] #Circular Code
setX

for (i in 1:length(codons)) {
  if (!partOfCircularCode(codons[[i]],setX)) { #not part of the Circular Code
    aminoAcid = translate(codons[[i]])
    codonsForAA = getCodesForAA(aminoAcid)
    
    ccCodonsForAA = getCircularCodes(codonsForAA, setX)
    if (!isEmpty(ccCodonsForAA)) { # Codons for a simple Substitution found
      newCodon = codonSubstitution(toString(codons[[i]]), ccCodonsForAA)
      
      codons[[i]] = newCodon #TODO: Codon ersetzten
      print(codons[[i]])
    } else { # no Codons found, amino acid has to be replaced
      print("kommt noch")
    }
  }
}

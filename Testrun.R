
library(Biostrings)
library(ccmotif)
source("ChangeCodons.R")

seqSet = readDNAStringSet("cds/ena-sars.fasta") # RNA sequence set
seq=seqSet[[11]] #Single RNA sequence
codonsOfSequence = codons(seq) # List of codons in rnaseq

setX = codes.c3[[26]] #Circular Code
setX

#subMatrix = BLOSUM62 #Substitution Matrix FUNKTIONIERT NOCH NICHT

newSeqString = "" #new sequence

for (i in 1:length(codonsOfSequence)) {
  newCodon = codonsOfSequence[[i]] #assume that no change is needed, new codon = old codon
  
  if (!partOfCircularCode(codonsOfSequence[[i]],setX)) { #not part of the Circular Code
    aminoAcid = translate(codonsOfSequence[[i]])
    codonsForAA = getCodesForAA(aminoAcid)
    ccCodonsForAA = getCircularCodes(codonsForAA, setX)
    
    if (!isEmpty(ccCodonsForAA)) { # there are other codons part of X coding for the same amino acid 
      newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsForAA) #replace codon
      #print(codonsOfSequence[[i]])
    } 
    else { # there are no other codons part of X coding for the same amino acid, try to replace amino acid
      newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, 0)

      if (newAminoAcid != AAString("")) { # there are >0 amino acids above threshold
        codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
        ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
        newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA)
      }
    }
  }
  newSeqString <<- paste(newSeqString, newCodon, sep="") #build new sequence
  newSeq = DNAString(newSeqString)
  
}

print(seq)
print("neu:")
print(newSeq)

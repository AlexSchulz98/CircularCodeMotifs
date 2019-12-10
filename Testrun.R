
library(Biostrings)
library(ccmotif)
source("ChangeCodons.R")
source("RawData.R")

seqSet = readDNAStringSet("cds/ena-sars.fasta") # RNA sequence set
seq=seqSet[[1]] #Single RNA sequence
codonsOfSequence = codons(seq) # List of codons in rnaseq

setX = codes.c3[[26]] #Circular Code

#subMatrix = BLOSUM62 #Substitution Matrix FUNKTIONIERT NOCH NICHT

newSeqString = "" #new sequence

for (i in 1:length(codonsOfSequence)) {
  newCodon = codonsOfSequence[[i]] #assume that no change is needed, new codon = old codon
  #print(paste("Ursprungscodon", newCodon))
  print(i)
  
  if (!partOfCircularCode(codonsOfSequence[[i]],setX)) { #not part of the Circular Code
    aminoAcid = translate(codonsOfSequence[[i]])
    codonsForAA = getCodesForAA(aminoAcid)
    ccCodonsForAA = getCircularCodes(codonsForAA, setX)
    
    #print("nicht Teil von X")
    
    if (!isEmpty(ccCodonsForAA)) { # there are other codons part of X coding for the same amino acid 
      newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsForAA) #replace codon
      #print(paste("andere Codons gefunden, das hier ist es geworden:", newCodon))
    } 
    else { # there are no other codons part of X coding for the same amino acid, try to replace amino acid
      newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, 0)
      
      #print("keine anderen Codons gefunden, Aminosäure wird ersetzt")

      if (newAminoAcid != AAString("")) { # there are >0 amino acids above threshold
        codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
        ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
        newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA)
        #print(paste("Ersetzung erfolgreich:", newCodon, "Coded für:", newAminoAcid, "(vorher", translate(codonsOfSequence[[i]]), ")"))
      }
    }
  }
  newSeqString <<- paste(newSeqString, newCodon, sep="") #build new sequence
  newSeq = DNAString(newSeqString)
  
}

# fileConn = file("output.txt") #newSequence to txt file
# writeLines(newSeq, fileConn) 
# close(fileConn)

print(setX)
print(toString(getAminoAcidsCodedByX(setX)))

print(seq)
print("neu:")
print(newSeq)

outputMatrix = codonCount(ar, seq,newSeq)

print(outputMatrix)


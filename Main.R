library(Biostrings)
library(ccmotif)
source("CodeManipulation.R")


#' generates modified sequence
#' @param seq single DNA sequence from sequence set, Format: seqSet[[i]]
#' @param setX circular code from codes.c3, Format: codes.c3[[i]]
#' @return new modified sequence
main = function(seq, setX) {
  
  codonsOfSequence = codons(seq) # List of codons in rnaseq
  #subMatrix = BLOSUM62 #Substitution Matrix 
  #TODO: BLOSUM selbst angeben FUNKTIONIERT NOCH NICHT
  
  newSeqString = "" #new sequence
  
  for (i in 1:length(codonsOfSequence)) {
    newCodon = codonsOfSequence[[i]] #assume that no change is needed, new codon = old codon
    #print(paste("Ursprungscodon", newCodon))
    #print(i)
    
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
        oldAminoAcid = translate(codonsOfSequence[[i]])
        newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, 0)
        
        #print("keine anderen Codons gefunden, Aminosäure wird ersetzt")
        
        if (newAminoAcid != oldAminoAcid) { # new amino acid found
          codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
          ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
          if (!isEmpty(ccCodonsFornewAA)) {
            newCodon <<- codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA)
            #print(paste("Ersetzung erfolgreich:", newCodon, "Coded für:", newAminoAcid, "(vorher", translate(codonsOfSequence[[i]]), ")"))
          }
        }
      }
    }
    newSeqString = paste(newSeqString, newCodon, sep="") #build new sequence
    #print(newSeqString)
  }
  
  newSeq = DNAString(newSeqString)
  
  return(newSeq)
  
}
# library(Biostrings)
# library(ccmotif)
#source("CodeManipulation.R")


#' generates modified sequence
#' @param seq single DNA sequence from sequence set, Format: seqSet[[i]]
#' @param setX circular code from codes.c3, Format: codes.c3[[i]]
#' @return new modified sequence
main = function(seq, setX) {

  codonsOfSequence = suppressWarnings(codons(seq))# List of codons in rnaseq
  #TODO: handle IUPAC codons
  newSeqString = ""
  
  for (i in 1:length(codonsOfSequence)) {
    newCodon = codonsOfSequence[[i]] #assume that no change is needed, new codon = old codon
    
    if (!partOfCircularCode(codonsOfSequence[[i]], setX)) {
      #not part of the Circular Code
      
      aminoAcid = Biostrings::translate(codonsOfSequence[[i]])
      codonsForAA = getCodesForAA(aminoAcid)
      ccCodonsForAA = getCircularCodes(codonsForAA, setX)
      
      if (!isEmpty(ccCodonsForAA)) {
        # there are other codons part of X coding for the same amino acid
        newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsForAA) #replace codon

      }
      else {
        # there are no other codons part of X coding for the same amino acid, try to replace amino acid
        
        oldAminoAcid = Biostrings::translate(codonsOfSequence[[i]])
        newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, 0)
        
        if (newAminoAcid != oldAminoAcid) {
          # new amino acid found
          codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
          ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
          if (!isEmpty(ccCodonsFornewAA)) {
            newCodon <<-
              codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA)
          }
        }
      }
    }

    newSeqString = paste(newSeqString, newCodon, sep = "") #build new sequence
  }
  
  newSeq = DNAString(newSeqString)
  
  return(newSeq)
}

library(Biostrings)
library(ccmotif)
source("ChangeCodons.R")
source("RawData.R")


#' writes new sequence into fasta file
#' @param seq single DNA sequence from sequence set, Format: seqSet[[i]]
#' @param setX circular code from codes.c3, Format: codes.c3[[i]]
#' @return new modified sequence
main = function(seq, setX) {

start_time <- Sys.time()

#seqSet = readDNAStringSet("cds/ena-sars.fasta") # RNA sequence set
#seq=seqSet[[11]] #Single RNA sequence
codonsOfSequence = codons(seq) # List of codons in rnaseq

#setX = codes.c3[[26]] #Circular Code

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
      newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, 0)
      
      #print("keine anderen Codons gefunden, Aminosäure wird ersetzt")

      if (newAminoAcid != AAString("")) { # there are >0 amino acids above threshold
        codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
        ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
        newCodon <<- codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA)
        #print(paste("Ersetzung erfolgreich:", newCodon, "Coded für:", newAminoAcid, "(vorher", translate(codonsOfSequence[[i]]), ")"))
      }
    }
  }
  newSeqString = paste(newSeqString, newCodon, sep="") #build new sequence
  #print(newSeqString)
}

newSeq = DNAString(newSeqString)

#outputMatrix = codonCount(ar, seq,newSeq)

#print(setX)
#print(toString(getAminoAcidsCodedByX(setX)))
#print(seq)
#print("neu:")
#print(newSeq)
#print(outputMatrix)

#write DNA Sequence into file
#a = append, w = new file


#write codon count into file
#saveRDS(object = outputMatrix, file = "outputMatrixSars.RDS")

end_time <- Sys.time()
print(end_time - start_time)

return(newSeq)

}

#---------------------------------------------


seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
setX = codes.c3[[26]] #Circular Code

ar = array(data = 0,
           dim = c(64,64,length(seqSet)),
           dimnames = list(CODONS,
                           CODONS,1:length(seqSet)))

start_time_global <- Sys.time()

for (z in 1:length(seqSet)) {
  
  newSeq= main(seqSet[[z]],setX)
  write.fasta(sequences = newSeq, names = paste("SARS VIRUS CDS NUMBER",z), file.out = "output_SarsVirus_CircularCode26.fasta", open = "a", nbchar = 60, as.string = FALSE)
  
  outputMatrix = codonCount(ar,length(seqSet), seqSet[[z]],newSeq, z)
  ar = outputMatrix
  
  print(paste("Done with set",z,"/",length(seqSet)))
}

saveRDS(object = outputMatrix, file = "outputMatrix_SarsVirus_CircularCode26.RDS")

end_time_global <- Sys.time()
print(end_time_global - start_time_global)








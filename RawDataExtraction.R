library(Biostrings)

#' generates 3D matrix with zeros
#' @param dim1 number
#' @param dim2 number
#' @param dim3 number of different sequences
#' @param dimnamesXY concat string names of axis x and y
#' @return 3darray
generateEmptyTable= function(dim1,dim2,dimnamesXY){
  
  ar = array(data = 0,
             dim = c(dim1,dim2),
             dimnames = list(dimnamesXY,
                             dimnamesXY))
return(ar)
}

#' fills table by counting how often codon A has been replaced with codon B
#' @param ar matrix 64x64
#' @param seqA DNAString from biological sequence
#' @param seqB DNAString from modified sequence
#' @return array of values
codonCount = function(ar, seqA, seqB){
  
  codonsA = codons(seqA)
  codonsB = codons(seqB)
  
  for (i in 1:length(codonsA)) {
    tmpArray =  array(data = 0,dim = c(64,64),dimnames = list(CODONS,CODONS))
    tmpArray[toString(codonsA[[i]]),toString(codonsB[[i]])] = 1
    ar = ar+tmpArray
  }
  return(ar)
}

#' fills table by counting how often codon A has been replaced with codon B
#' @param ar matrix 21x21
#' @param seqA DNAString from biological sequence
#' @param seqB DNAString from modified sequence
#' @return array of values
aminoCount = function(ar, seqA, seqB){
  
  aminoA = translate(seqA)
  aminoB = translate(seqB)
  
  for (i in 1:length(aminoA)) {
    
    tmpArray =  array(data = 0,dim = c(21,21),dimnames = list(AMINOACIDS,AMINOACIDS))
    tmpArray[toString(aminoA[i]),toString(aminoB[i])] = 1
    ar = ar+tmpArray
  }
  return(ar)
}

#' Amount of changed bases in a sequence compared to another
#' @param seqA original sequence
#' @param seqB modified sequence
baseCount = function(seqA, seqB){
  
  sum = 0
  
  for (i in 1:length(seqA)) {
    if (seqA[i] != seqB[i]) {
      sum = sum +1
    }
  }
  return(sum)
}

#' creates a binary sequence depending on motif (=1) and non-motif(=0) codons
#' @param seq dna sequence
#' @param setX circular codes
#' @return array cosisting of zeros and ones
sequenceToBinary = function(seq, setX){
  
  binSeq = array()
  
  codon = codons(seq)
  for (i in 1:length(codon)) {
    if (partOfCircularCode(codon[[i]],setX)) {
      binSeq[i] = 1
    } else {
      binSeq[i] = 0
    }
  }
  return(binSeq)
}

#'
#' @param binSeq array 
# plotBinarySequence = function(binSeq){
#   
#   motifs = factor(binSeq)
#   print(length(motifs))
#   position_in_sequence = length(binSeq)
#   print(position_in_sequence)
#   
#   cdplot(motifs ~ position_in_sequence) #TODO: Error
# }








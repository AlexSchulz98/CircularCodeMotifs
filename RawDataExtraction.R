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



#' Finds codons which could not be changed although they were not part of the circular code
#' @param table codon table
#' @param codes circular code
#' @return amout of unchanged non-circular code codons
unchangednonCCCodons = function(table, codes) {
  
  sum = 0
  dia = diag(table)
  
  for (i in 1:length(dia)) {
    label = names(dia[i])
    bool = partOfCircularCode(label, codes)
    
    if (bool == FALSE) {
      sum = sum + dia[[i]]
    }
  }
  
  return(sum)
}


#' assess how much a sequence has been changed by giving mutations scores
#' +1 for point mutations
#' +3 for changes of 1 or more basis but maintaining the amino acid
#' +5 for the change of an amino acid 
#' @param table equal sided table
#' @param codes circular code
valueTable = function(table, codes){
  
  sum = 0
  
  unc = unchangednonCCCodons(table, codes)
  sum = -0.5*unc
  
  dimensions = dim(table)
  dimNames = dimnames(table)
  dimNames1 = dimNames[[1]]
  dimNames2 = dimNames[[2]]
  
  for (i in 1:dimensions[1]) {
    
    codon1 = dimNames1[i]
    
    for (j in 1:dimensions[2]) {
      
      codon2 = dimNames2[j]
      
      if (codon1 != codon2) { #ignore diagonal
        
        tmp = compareCodons(codon1,codon2)
        sum = sum + tmp
      }
    }
  }
  
  return(sum)
}


#' compares two codons and scores them
#' Base 1 the same +4
#' Base 2 the same +4
#' Base 3 the same +2
#' 1/2/3 not the same, but purin/purin or pyrimidin/pyrimidin change +1 (amineChange function)
compareCodons = function(codon1, codon2) {
  
  split1 = unlist(strsplit(codon1, ""))
  split2 = unlist(strsplit(codon2, ""))
  
  score = 0
  
  if (split1[1] == split2[1]) {
    score = score + 4
  } else {
    amineCheckScore = amineChange(split1[1], split2[1])
    score = score + amineCheckScore
  }
  if (split1[2] == split2[2]) {
    score = score + 4
  } else {
    amineCheckScore = amineChange(split1[2], split2[2])
    score = score + amineCheckScore
  }
  if (split1[3] == split2[3]) {
    score = score + 2
  } else {
    amineCheckScore = amineChange(split1[3], split2[3])
    score = score + amineCheckScore
  }
  return(score)
}



#for (i in 1:HIERFEHLTWAS) {
  #Änderung aus tabelle berechenen
  #alles außerhalb der diagnole wert zuteilen
  # je nach änderung
   #aminostausch erkennung einbauen?
#}

#' Codon teil von X --> lassen?
#' Codons nicht mutierbar --> -1 ?
#' codons nach schema f bewerten (max 10 wird nie erreicht, 9 kann aber vorkommen, mittelwert 5?)


















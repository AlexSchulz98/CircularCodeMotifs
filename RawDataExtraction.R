library(Biostrings)

#' generates 3D matrix with zeros
#' @param dim1 number
#' @param dim2 number
#' @param dim3 number of different sequences
#' @param dimnamesXY string names of axis x and y
#' @return 3darray
generateEmptyTable= function(dim1,dim2,dim3,dimnamesXY){
  
  ar = array(data = 0,
             dim = c(dim1,dim2,dim3),
             dimnames = list(dimnamesXY,
                             dimnamesXY,1:dim3))
return(ar)
}

#' fills table by counting how often codon A has been replaced with codon B
#' @param ar matrix
#' @param dim3 3D depth of array
#' @param codonA DNAString from biological sequence
#' @param codonB DNAString from modified sequence
#' @param z current depth
codonCount = function(ar, dim3, seqA, seqB, z){
  
  codonsA = codons(seqA)
  codonsB = codons(seqB)
  
  for (i in 1:length(codonsA)) {
    tmpArray =  array(data = 0,dim = c(64,64,dim3),dimnames = list(CODONS,CODONS,1:dim3))
    tmpArray[toString(codonsA[[i]]),toString(codonsB[[i]]),z] = 1
    ar = ar+tmpArray
  }
  return(ar)
}

#' fills table by counting how often codon A has been replaced with codon B
#' @param ar matrix
#' @param dim3 3D depth of array
#' @param codonA DNAString from biological sequence
#' @param codonB DNAString from modified sequence
#' @param z current depth
aminoCount = function(ar, dim3, seqA, seqB, z){
  
  aminoA = translate(seqA)
  aminoB = translate(seqB)
  
  for (i in 1:length(aminoA)) {
    print(i)
    
    tmpArray =  array(data = 0,dim = c(21,21,dim3),dimnames = list(AMINOACIDS,AMINOACIDS,1:dim3))
    tmpArray[toString(aminoA[i]),toString(aminoB[i]),z] = 1
    ar = ar+tmpArray
  }
  return(ar)
}

#' Anzahl einzelener Mutationen
baseCount = function(codonB){
  
}

#' Durch Aminosäuren nach Anzahl Austauschen suchen
aminoAcidCount = function(codonA, codonB){
  
}

#' Motiflängen mit 1 und 0
motifLenght = function(codonA, codonB){
  
}



#' Anzahl ausgetauschter Codons --> Matrix?
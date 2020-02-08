library(Biostrings)
source("CodeManipulation.R")
source("Parameter.R")

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
  
  for (i in 1:length(codonsB)) {
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
  
  aminoA = Biostrings::translate(seqA)
  aminoB = Biostrings::translate(seqB)
  
  for (i in 1:length(aminoB)) {
    
    tmpArray =  array(data = 0,dim = c(21,21),dimnames = list(AMINOACIDS,AMINOACIDS))
    tmpArray[toString(aminoA[i]),toString(aminoB[i])] = 1
    ar = ar+tmpArray
  }
  return(ar)
}



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


#' asseses how much a sequence has been changed by giving mutations scores
#' Rating for the modified sequences:
#' Base 1 the same +4
#' Base 2 the same +4
#' Base 3 the same +2
#' 1/2/3 not the same, but purin/purin or pyrimidin/pyrimidin change +1
#' Sum of this rating is divided by the number of codons to get the average
#' because of the rating range from 0 to 10 the score is divided by 10 to get a normalized scale from 0 to 1
#' @param table codon changes
#' @return score between 0 and 1. 1=good/minor mutations. 
getEditScore = function(table){
  
  score = 0 
  
  total = sum(table) # codons in this sequence
  unc = sum(diag(table)) # unchanged codons
  sumCC = total-unc # changed codons
  

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
        tmp = tmp*table[codon1,codon2] # multiply score with amount of changes (which are stored in table)
        score = score + tmp
      }
    }
  }
  
  result = score/sumCC # calc average (score is between 0-10)
  result = result/10 # normalize to 0-1
  
  return(result)
}

#' determines the distance of a sequence to a perfect (full circular code) sequence
#' @param table codon table
#' @param editScore calculated edit score between 0 and 1
#' @param code circular code used
#' @return score between 0 and 1. 0=good/close to perfect sequence.
getEditDistance = function(table, editScore, code){
  
  total = sum(table) # codons in this sequence
  unc = sum(diag(table)) # unchanged codons
  
  cc = total-unc # changed codons
  ncp = unchangednonCCCodons(table,code) # no change possible
  
  #sumCCandNCP = sumCC+ncp # non cc codons in the original sequence
  
  #cc_p = cc/sumCCandNCP #percent that has been changed
  #ncp_p = ncp/sumCCandNCP #percent that could not be changed
  cc_p = cc/total # percent that has been changed
  ncp_p = ncp/total #percent that could not be changed
  
  distance = cc_p*(1-editScore) + ncp_p 
  
  return(distance)
  
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


















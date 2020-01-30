library(Biostrings)

# 64 codons
CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
# 20 amino acids + 1 *  for stop codons
AMINOACIDS = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

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


#' assess how much a sequence has been changed by giving mutations scores
#' Rating for codons that did not need to be changed, because they already were part of the circular code: 10
#' Rating for the modified sequences:
#' Base 1 the same +4
#' Base 2 the same +4
#' Base 3 the same +2
#' 1/2/3 not the same, but purin/purin or pyrimidin/pyrimidin change +1
#' Sum of both ratings is multiplied with the amount of circular codes in the new sequence (code usage). 
#' Logarithm to center score around 0. Max. score is 1. Worst score is -1
#' @param table equal sided table
#' @param codes circular code
#' @param cu code usage
#' @return score between -1 and 1. Best is 1. 
valueTable = function(table, codes, cu){
  
  score = 0
  total = sum(table)
  dia = sum(diag(table))
  
  unc = unchangednonCCCodons(table, codes)
  
  sumCC = total-unc
  ch = dia-unc
  score = ch*10
  

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
        score = score + tmp
      }
    }
  }
  
  result = score/sumCC
  result = result*cu
  result = log(result,base = 10)
  
  return(result)
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


















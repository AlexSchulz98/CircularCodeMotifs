library(Biostrings)
source('D:/Projekte/BA Circular Code/Scripts/Sequences.R', echo=TRUE)
source('D:/Projekte/BA Circular Code/Scripts/Mutations.R', echo=TRUE)

# To name rows and columns of matrix:
# 64 codons
CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
# 20 amino acids + 1 *  for stop codons
AMINOACIDS = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")


#' Generate 3D matrix with zeros
#'
#' @param dim1 number
#' @param dim2 number
#' @param dimnamesXY vector String names
#'
#' @return array
#' @export
#'
#' @examples
#' generateEmptyTable(64, 64, CODONS)
generateEmptyTable= function(dim1,dim2,dimnamesXY){
  
  ar = array(data = 0,
             dim = c(dim1,dim2),
             dimnames = list(dimnamesXY,
                             dimnamesXY))
return(ar)
}


#' Fill codon matrix
#'
#' @param ar matrix 64x64
#' @param seqA DNAString sequence
#' @param seqB DNAString sequence
#'
#' @return matrix 64x64
#' @export
#'
#' @examples
#' ar = generateEmptyTable(64, 64, CODONS)
#' seqA = DNAString("ATGATGATG")
#' seqB = DNAString("ATGATGATG")
#' codonCount(ar, seqA, seqB]
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

#' Fill amino matrix (deprecated)
#'
#' @param ar matrix 21x21
#' @param seqA DNAString sequence
#' @param seqB DNAString sequence
#'
#' @return matrix 21x21
#' @export
#'
#' @examples
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



#' Amount of codons change not possible (on diagonal but not part of x)
#'
#' @param table matrix 
#' @param codes codes.code circular code
#'
#' @return number
#' @export
#'
#' @examples
#' unchangednonCCCodons(rds,codes.c3[[23]])
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


#' Get Edit-Score by 4-4-2-1 Scheme
#'
#' @param table (filled) codon matrix 64x64
#'
#' @return score between 0 and 1
#' @export
#'
#' @examples
#' getEditScore(rds)
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


#' Get Edit_Distance
#'
#' @param table matrix 64x64
#' @param editScore  edit score between 0 and 1
#' @param code codes.code circular code
#'
#' @return score between 0 and 1
#' @export
#'
#' @examples
#' getEditDistance(rds,editScore,codes.c3[[23]])
getEditDistance = function(table, editScore, code){
  
  total = sum(table) # codons in this sequence
  unc = sum(diag(table)) # unchanged codons
  
  cc = total-unc # changed codons
  ncp = unchangednonCCCodons(table,code) # no change possible
  
  cc_p = cc/total # percent that has been changed
  ncp_p = ncp/total #percent that could not be changed
  
  distance = cc_p*(1-editScore) + ncp_p 
  
  return(distance)
  
}


#' Simple comparision with 4-4-2-1 scheme
#'
#' @param codon1 String Codon
#' @param codon2 String Codon
#'
#' @return score (0-10)
#' @export
#'
#' @examples
#' compareCodons("AAA","AAT")
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



















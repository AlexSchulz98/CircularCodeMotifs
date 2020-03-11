library(Biostrings)
source('../Scripts/Mutations.R', echo=TRUE)

#' Change Reading Frame (deletes first 1 or 2 bases)
#'
#' @param targetFrame
#' @param seq DNAString sequence
#'
#' @return DNAString sequence
#' @export
#'
#' @examples
#' seq = DNAString("ATGCTGTTC")
#' targetFrame = 1
#' changeReadingFrame(targetFrame, seq)
changeReadingFrame = function(targetFrame, seq) {
  if (targetFrame < 0)
    stop("target frame cannot be lower than 0")
  
  return(subseq(seq, start = targetFrame + 1))
}


#' Get amount of CDS with IUPAC symbols (except A-G-C-T)
#'
#' @param seqSet DNAStringSet
#' @param rangeStart number, CDS to start with
#' @param rangeEnd number, CDS to end with
#'
#' @return number
#' @export
#'
#' @examples
#' seqSet = readDNAStringSet("../BA Circular Code/cds/ena-sars.fasta")
#' rangeStart = 1
#' rangeEnd = 11
#' getIUPACCodesCount(seqSet, rangeStart, rangeEnd)
getIUPACCodesCount = function(seqSet, rangeStart = 1, rangeEnd) {
  sum = 0
  
  for (i in rangeStart:rangeEnd) {
    skip = FALSE
    tryCatch({
      tmp = suppressWarnings(Biostrings::codons(seqSet[[i]]))
    }, error = function(e) {
      skip <<- TRUE
      sum <<- sum + 1
    })
    if (skip == TRUE) {
      next
    }
  }
  return(sum)
}


#' delete sequences with IUPAC codes
#'
#' @param seqSet DNAStringSet
#'
#' @return DNAStringSet
#' @export
#'
#' @examples
#' seqSet = readDNAStringSet("../BA Circular Code/cds/ena-sars.fasta")
#' deleteIUPACSequences(seqSet)
deleteIUPACSequences = function(seqSet) {
  newSeqSet = DNAStringSet()
  j = 1
  
  for (i in 1:length(seqSet)) {
    skip = FALSE
    tryCatch({
      tmp = suppressWarnings(Biostrings::codons(seqSet[[i]]))
    }, error = function(e) {
      skip <<- TRUE
    })
    if (skip == TRUE) {
      next
    }
    else {
      newSeqSet[j] = seqSet[i]
      j = j + 1
    }
    
  }
  return(newSeqSet)
}

#' Shift each codon of a code
#'
#' @param code number
#' @param shift 1 or 2
#'
#' @return vector with String codons
#' @export
#'
#' @examples
#' code = 23
#' shift = 1
#' shiftCode(code,shift)
shiftCode = function(code, shift) {
  mainCode = codes.c3[[code]]
  mainCode = mainCode[[2]]
  
  newCode = c()
  
  
  for (i in 1:length(mainCode)) {
    newCodon = ""
    
    codon = mainCode[i]
    
    codon1 = unlist(strsplit(codon, ""))
    
    if (shift == 1) {
      newCodon = paste(codon1[2], codon1[3], codon1[1], sep = "")
    } else {
      newCodon = paste(codon1[3], codon1[1], codon1[2], sep = "")
    }
    
    newCode[i] = newCodon
    
  }
  
  return(newCode)
  
}


#' Check vector of codons if in code
#'
#' @param codons vector of String codons
#' @param code number
#'
#' @return vector with FALSE / TRUE
#' @export
#'
#' @examples
#' codons = c("ATG","ATT","AGG")
#' code = 23
#' partOfCircularCode_List(codons, code)
partOfCircularCode_List = function(codons, code) {
  result = c()
  
  
  for (i in 1:length(codons)) {
    tmp = partOfCircularCode(codons[i], codes.c3[[code]])
    
    result[i] = tmp
  }
  
  return(result)
  
}




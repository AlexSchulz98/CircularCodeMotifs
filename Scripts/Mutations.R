library(Biostrings)
library(ccmotif)

#' Check if codon is part of the circular code
#'
#' @param codon 3-letter DNAString codon
#' @param code codes.code circular code
#'
#' @return TRUE / FALSE
#' @export
#'
#' @examples 
#' codon = DNAString("ATT")
#' code = codes.c3[[23]]
#' partOfCircularCode(codon,code)
partOfCircularCode = function(codon, code){
  
  if (grepl(codon, code[2])) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

#' Codons for amino acids
#'
#' @param aminoAcid AAString amino acid
#'
#' @return vector with String codons
#' @export
#'
#' @examples
#' aminoAcid = AAString("A")
#' getCodesForAA(aminoAcid)
getCodesForAA = function(aminoAcid) {
  if (aminoAcid == AAString("A")) {
    return(c("GCC", "GCT", "GCG", "GCA"))
  } else if (aminoAcid == AAString("R")) {
    return(c("CGA", "CGG", "CGT", "CGC", "AGG", "AGA"))
  } else if (aminoAcid == AAString("N")) {
    return(c("AAT", "AAC"))
  } else if (aminoAcid == AAString("D")) {
    return(c("GAT", "GAC"))
  } else if (aminoAcid == AAString("C")) {
    return(c("TGT", "TGC"))
  } else if (aminoAcid == AAString("Q")) {
    return(c("CAA", "CAG"))
  } else if (aminoAcid == AAString("E")) {
    return(c("GAA", "GAG"))
  } else if (aminoAcid == AAString("G")) {
    return(c("GGA", "GGG", "GGT", "GGC"))
  } else if (aminoAcid == AAString("H")) {
    return(c("CAT", "CAC"))
  } else if (aminoAcid == AAString("I")) {
    return(c("ATA", "ATT", "ATC"))
  } else if (aminoAcid == AAString("L")) {
    return(c("TTA", "TTG", "CTA", "CTG", "CTT", "CTC"))
  } else if (aminoAcid == AAString("K")) {
    return(c("AAA", "AAG"))
  } else if (aminoAcid == AAString("M")) {
    return(c("ATG"))
  } else if (aminoAcid == AAString("F")) {
    return(c("TTT", "TTC"))
  } else if (aminoAcid == AAString("P")) {
    return(c("CCA", "CCT", "CCG", "CCC"))
  } else if (aminoAcid == AAString("S")) {
    return(c("AGC", "AGT", "TCC", "TCT", "TCG", "TCA"))
  } else if (aminoAcid == AAString("T")) {
    return(c("ACC", "ACA", "ACG", "ACT"))
  } else if (aminoAcid == AAString("W")) {
    return(c("TGG"))
  } else if (aminoAcid == AAString("Y")) {
    return(c("TAC", "TAT"))
  } else if (aminoAcid == AAString("V")) {
    return(c("GTA", "GTT", "GTG", "GTC"))
  } else if (aminoAcid == AAString("*")) {
    return(c("TAA", "TAG", "TGA"))
  }
}

#' Filter non-cc codons
#'
#' @param codons vector of String codons
#' @param setX codes.code circular code
#'
#' @return vector with String codons
#' @export
#'
#' @examples
#' codons = c("GCC","GCT","GCG","GCA","ATT")
#' setX = code = codes.c3[[23]]
#' getCircularCodes(codons, setX)
getCircularCodes = function(codons, setX){
  x = c() 
  for (i in 1:length(codons)) {
    x[i] = partOfCircularCode(codons[i], setX)
  }
  return(codons[x])
}


#' Purin/Pyrimidin part of the 4-4-2-1 scheme
#' A <--> G || C <--> T
#'
#' @param oldBase String single base
#' @param newBase String single base
#'
#' @return 0 or 1
#' @export
#'
#' @examples
#' oldBase = "A"
#' newBase = "G"
#' amineChange(oldBase, newBase)
amineChange = function(oldBase, newBase){
  if ((oldBase == "G"|| oldBase == "A") && (newBase == "G" || newBase == "A")) {
    return(1) #pyrimidin replaced with pyrimidin
  } else if((oldBase == "C"|| oldBase == "T") && (newBase == "C" || newBase == "T")) {
    return(1)
  } else {
    return(0)
  }
}


#' Select best codon with 4-4-2-1 scheme
#'
#' @param oldCodon String codon
#' @param newCodons vector with String codons
#'
#' @return String codon
#' @export
#'
#' @examples
#' oldCodon = "AAA"
#' newCodons = c("AAG","AAT","AGA","GAA")
#' codonSubstitution(oldCodon, newCodons)
codonSubstitution = function(oldCodon, newCodons) {
  score = -1
  oldBases = unlist(strsplit(oldCodon, ""))
  for (i in 1:length(newCodons)) {
    tmpscore = 0
    newBases = unlist(strsplit(newCodons[i], ""))
    if (oldBases[1] == newBases[1]) {
      tmpscore = tmpscore + 4
    } else {
      amineCheckScore = amineChange(oldBases[1], newBases[1])
      tmpscore = tmpscore + amineCheckScore
    }
    if (oldBases[2] == newBases[2]) {
      tmpscore = tmpscore + 4
    } else {
      amineCheckScore = amineChange(oldBases[2], newBases[2])
      tmpscore = tmpscore + amineCheckScore
    }
    if (oldBases[3] == newBases[3]) {
      tmpscore = tmpscore + 2
    } else {
      amineCheckScore = amineChange(oldBases[3], newBases[3])
      tmpscore = tmpscore + amineCheckScore
    }
    
    if (tmpscore > score) {
      score = tmpscore
      newCodon = newCodons[i]
    }
  }
  
  return(newCodon)
}

#' Amino acids coded by a vector of codons
#'
#' @param setX codes.code circular code
#'
#' @return AAStringSet with (unique) amino acids
#' @export
#'
#' @examples
#' setX = codes.c3[[23]]
#' getAminoAcidsCodedByX(setX)
getAminoAcidsCodedByX = function(setX) {
  xCodes = strsplit(setX[[2]], split = " ")
  xCodeAminoAcids = vector("list", length(xCodes))
  
  for (i in 1:length(xCodes)) {
    xCodeDNA = DNAString(xCodes[[i]])
    xCodeAminoAcids[[i]] = Biostrings::translate(xCodeDNA)
  }
  xCodeAminoAcids = unique(AAStringSet(xCodeAminoAcids))
  return(xCodeAminoAcids)
}


#' Select best codon with BLOSUM62 matrix
#'
#' @param codon 3-letter DNAString codon
#' @param setX codes.code circular code
#' @param threshold number of lowest possible score for eligible substitution
#'
#' @return AAString amino acid
#' @export
#'
#' @examples
#' codon = DNAString("AAA")
#' setX = codes.c3[[23]]
#' threshold = 0
#' aminoAcidSubstitution(codon, setX, threshold)
aminoAcidSubstitution = function(codon, setX, threshold = 0) {
  data("BLOSUM62")
  xCodeAminoAcids = getAminoAcidsCodedByX(setX)
  
  aminoAcid = Biostrings::translate(codon)
  bestAminoAcid = aminoAcid # old aa in case no match is found
  
  for (i in 1:length(xCodeAminoAcids)) {
    score = BLOSUM62[toString(aminoAcid), toString(xCodeAminoAcids[[i]])]
    
    if (score >= threshold) {
      bestAminoAcid = xCodeAminoAcids[[i]]
      threshold = score
    }
  }
  return(bestAminoAcid)
}

#' Mutate a sequence with a circular code
#'
#' @param seq DNAString sequence
#' @param setX codes.code circular code
#'
#' @return DNAString sequence
#' @export
#'
#' @examples
#' seqSet = readDNAStringSet("../BA Circular Code/cds/ena-sars.fasta")
#' seq = seqSet[[11]]
#' setX = codes.c3[[23]]
#' subAminos = TRUE
#' threshold = 0
#' mutateSequence(seq,setX, threshold, subAminos)
mutateSequence = function(seq, setX, subAminos = FALSE, threshold = 0 ) {
  codonsOfSequence = suppressWarnings(codons(seq)) # warning if seq not a multiple of 3
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
        if (subAminos == TRUE) {
          # there are no other codons part of X coding for the same amino acid, try to replace amino acid
          
          oldAminoAcid = Biostrings::translate(codonsOfSequence[[i]])
          newAminoAcid = aminoAcidSubstitution(codonsOfSequence[[i]], setX, threshold)
          
          if (newAminoAcid != oldAminoAcid) {
            codonsFornewAA = getCodesForAA(newAminoAcid) #get codons for this amino acid
            ccCodonsFornewAA = getCircularCodes(codonsFornewAA, setX) #get cc codons for this amino acid
            if (!isEmpty(ccCodonsFornewAA)) {
              newCodon = codonSubstitution(toString(codonsOfSequence[[i]]), ccCodonsFornewAA) 
            }
          }
        }
      }
    }
    
    newSeqString = paste(newSeqString, newCodon, sep = "") #build new sequence
  }
  
  newSeq = DNAString(newSeqString)
  
  return(newSeq)
}


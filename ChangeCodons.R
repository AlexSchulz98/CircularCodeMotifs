library(Biostrings)

#' Check if codon is part of the circular codes
#' @param codon xStringViews or String codon, single
#' @param codes format codes.c3[[i]]) circular code
#' codes[1]: code number
#' codes[2]: List of Codons
partOfCircularCode = function(codon, codes){
  
  if (grepl(codon, codes[2])) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

#' replace a single codon
#' @param oldCodon String codon
#' @param newCodons concat with String entries
#' change purin with purin and pyrimidin with pyrimidin
#' use as few mutations as possible
#' mutations at the 3rd base are scored better (Wobble Hypothesis)
#' TODO: Score überarbeiten
#' @return String codon
codonSubstitution = function(oldCodon, newCodons){
  score = -1
  oldBases = unlist(strsplit(oldCodon,""))
  for (i in 1:length(newCodons)) {
    tmpscore = 0
    newBases = unlist(strsplit(newCodons[i],""))
    if (oldBases[1] == newBases[1]) { 
      tmpscore = tmpscore+3
    }
    if (oldBases[2] == newBases[2]) { 
      tmpscore = tmpscore+3
    }
    if (oldBases[3]==newBases[3]) {
      tmpscore = tmpscore+2
    }
    #change Purin with Purin and Pyrimidin with Pyrimidin 
    if((oldBases[3] == "G"|| oldBases[3] == "A") && (newBases[3] == "G" || newBases[3] == "A")){
      tmpscore = tmpscore+1 #pyrimidin replaced with pyrimidin
    } else if ((oldBases[3] == "C"|| oldBases[3] == "T") && (newBases[3] == "C" || newBases[3] == "T")) {
      tmpscore = tmpscore+1 #purin replaced with purin
    }
    if (tmpscore > score) {
      score = tmpscore
      newCodon = newCodons[i]
    }
  }

  return(newCodon)
}

#' Replace amino acid based on substitution matrix
#' TODO: @param subMatrix matrix. Datenbank wird nicht gefunden, wenn man sie einer variablen übergibt
#' @param codon codon in old sequence
#' @param setX circular code
#' @param threshold lowest possible score in matrix for eligible replacement
#' @return AAString amino acid
aminoAcidSubstitution = function(codon, setX, threshold){
  
  data("BLOSUM62")
  xCodeAminoAcids = getAminoAcidsCodedByX(setX)
  
  aminoAcid = translate(codon)
  
  newAminoAcid = AAString("")
  
  for (i in 1:length(xCodeAminoAcids)) {
    score = BLOSUM62[toString(aminoAcid),toString(xCodeAminoAcids[[i]])]
    
    if (score >= threshold) {
      newAminoAcid = xCodeAminoAcids[[i]]
      threshold=score
    }
  }
  return(newAminoAcid)
}

#' returns a list of amino acids which are coded by codons part of X
#' @param setX circular code
#' @return AAStringSet, unique
getAminoAcidsCodedByX = function(setX) {
  xCodes = strsplit(setX[[2]], split=" ")
  xCodeAminoAcids = vector("list", length(xCodes))
  
  for (i in 1:length(xCodes)) {
    xCodeDNA = DNAString(xCodes[[i]])
    xCodeAminoAcids[[i]] = translate(xCodeDNA)
  }
  xCodeAminoAcids = unique(AAStringSet(xCodeAminoAcids))
  return(xCodeAminoAcids)
}

#' takes a list of cc- and non cc-codons and return s only those which are part of the Circular Code
#' @param codons list of cc- and non cc-codons
#' @return list of circular code Codons
getCircularCodes = function(codons, setX){
  x = c() #boolean vector to delete non-cc codes from codons
  for (i in 1:length(codons)) {
    x[i] = partOfCircularCode(codons[i], setX)
  }
  return(codons[x])
}

#' all possible codons fpr each amino acid or stop codon
#' @param aminoAcid AAString amino acid
#' @return concat with String entries
getCodesForAA = function(aminoAcid){
  if (aminoAcid == AAString("A")) {
    return(c("GCC","GCT","GCG","GCA"))
  } else if (aminoAcid == AAString("R")) {
    return(c("CGA","CGG","CGT","CGC","AGG","AGA"))
  } else if (aminoAcid == AAString("N")) {
    return(c("AAT","AAC"))
  } else if (aminoAcid == AAString("D")) {
    return(c("GAT","GAC"))
  } else if (aminoAcid == AAString("C")) {
    return(c("TGT","TGC"))
  } else if (aminoAcid == AAString("Q")) {
    return(c("CAA","CAG") )
  } else if (aminoAcid == AAString("E")) {
    return(c("GAA","GAG"))
  } else if (aminoAcid == AAString("G")) {
    return(c("GGA","GGG","GGT","GGC"))
  } else if (aminoAcid == AAString("H")) {
    return(c("CAT","CAC"))
  } else if (aminoAcid == AAString("I")) {
    return(c("ATA","ATT","ATC"))
  } else if (aminoAcid == AAString("L")) {
    return(c("TTA","TTG","CTA","CTG","CTT","CTC"))
  } else if (aminoAcid == AAString("K")) {
    return(c("AAA","AAG"))
  } else if (aminoAcid == AAString("M")) {
    return(c("ATG"))
  } else if (aminoAcid == AAString("F")) {
    return(c("TTT","TTC"))
  } else if (aminoAcid == AAString("P")) {
    return(c("CCA","CCT","CCG","CCC"))
  } else if (aminoAcid == AAString("S")) {
    return(c("AGC","AGT","TCC","TCT","TCG","TCA"))
  } else if (aminoAcid == AAString("T")) {
    return(c("ACC","ACA","ACG","ACT"))
  } else if (aminoAcid == AAString("W")) {
    return(c("TGG"))
  } else if (aminoAcid == AAString("Y")) {
    return(c("TAC","TAT"))
  } else if (aminoAcid == AAString("V")) {
    return(c("GTA","GTT","GTG","GTC"))
  } else if (aminoAcid == AAString("*")) {
    return(c("TAA","TAG","TGA"))
  }
}
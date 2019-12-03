library(Biostrings)

#' Was ich habe:
#' Codes.C3[[i]] --> Liste mit 2 elementen: in [1] Beschreibung und in [2] die Codes
#' seq --> DNAString


#' Check if codon is part of the circular codes
#' @param codon single codon as xStringViews or as String
#' @param codes Circular code as single list (Format: codes.c3[[i]])
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

#' Dient dazu ein Codon auszutauschen
#' @param test nutzlos
#' Beachten: Purin vs Pyridin
#' Möglichst wenig mutationen

codonSubstitution = function(codon){
  
  
  return(newCodon)
}


#' Dient dazu eine einzlene Base auszutauschen

baseSubstitution = function(kommtmorgen){
  
}





getCodesForAA = function(){
  
}


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
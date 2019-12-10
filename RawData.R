library(Biostrings)

CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")


ar = array(data = 0,
           dim = c(64,64),
           dimnames = list(CODONS,
                           CODONS))

#' fills table by counting how often codon A has been replaced with codon B
#' @param codonA DNAString from biological sequence
#' @param codonB DNAString from modified sequence
codonCount = function(ar, seqA, seqB){
  codonsA = codons(seqA)
  codonsB = codons(seqB)
  print(length(codonsA))
  
  for (i in 1:length(codonsA)) {
    tmpArray =  array(data = 0,dim = c(64,64),dimnames = list(CODONS,CODONS))
    print(toString(codonsA[[i]]))
    print(toString(codonsB[[i]]))
    tmpArray[toString(codonsA[[i]]),toString(codonsB[[i]])] = 1
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
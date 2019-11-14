library(Biostrings)

#amino acids with all codes
ala = c("GCC","GCT","GCG","GCA")
arg = c("CGA","CGG","CGT","CGC","AGG","AGA")
asn = c("AAT","AAC")
asp = c("GAT","GAC")
cys = c("TGT","TGC")
gln = c("CAA","CAG") 
glu = c("GAA","GAG")
gly = c("GGA","GGG","GGT","GGC")
his = c("CAT","CAC")
ile = c("ATA","ATT","ATC")
leu = c("TTA","TTG","CTA","CTG","CTT","CTC")
lys = c("AAA","AAG")
met = c("ATG")
phe = c("TTT","TTC")
pro = c("CCA","CCT","CCG","CCC")
ser = c("AGC","AGT","TCC","TCT","TCG","TCA")
thr = c("ACC","ACA","ACG","ACT")
trp = c("TGG")
tyr = c("TAC","TAT")
val = c("GTA","GTT","GTG","GTC")

getCodesForAA = function(aminoAcid){
  if (aminoAcid == "A") {
    return(ala)
  } else if (aminoAcid == "R") {
    return(arg)
  } else if (aminoAcid == "N") {
    return(asn)
  } else if (aminoAcid == "D") {
    return(asp)
  } else if (aminoAcid == "C") {
    return(cys)
  } else if (aminoAcid == "Q") {
    return(gln)
  } else if (aminoAcid == "E") {
    return(glu)
  } else if (aminoAcid == "G") {
    return(gly)
  } else if (aminoAcid == "H") {
    return(his)
  } else if (aminoAcid == "I") {
    return(ile)
  } else if (aminoAcid == "L") {
    return(leu)
  } else if (aminoAcid == "K") {
    return(lys)
  } else if (aminoAcid == "M") {
    return(met)
  } else if (aminoAcid == "F") {
    return(phe)
  } else if (aminoAcid == "P") {
    return(pro)
  } else if (aminoAcid == "S") {
    return(ser)
  } else if (aminoAcid == "T") {
    return(thr)
  } else if (aminoAcid == "W") {
    return(trp)
  } else if (aminoAcid == "Y") {
    return(tyr)
  } else if (aminoAcid == "V") {
    return(val)
  }
}


getCodesOfInterest = function(circularCode){
  
  aminoAcidSeq = translate(circularCode)
  print(aminoAcidSeq)
  
}







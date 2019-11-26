library(Biostrings)
library(ccmotif)

#----
#amino acids with all codes
#ala = c("GCC","GCT","GCG","GCA")
#arg = c("CGA","CGG","CGT","CGC","AGG","AGA")
#asn = c("AAT","AAC")
#asp = c("GAT","GAC")
#cys = c("TGT","TGC")
#gln = c("CAA","CAG") 
#glu = c("GAA","GAG")
#gly = c("GGA","GGG","GGT","GGC")
#his = c("CAT","CAC")
#ile = c("ATA","ATT","ATC")
#leu = c("TTA","TTG","CTA","CTG","CTT","CTC")
#lys = c("AAA","AAG")
#met = c("ATG")
#phe = c("TTT","TTC")
#pro = c("CCA","CCT","CCG","CCC")
#ser = c("AGC","AGT","TCC","TCT","TCG","TCA")
#thr = c("ACC","ACA","ACG","ACT")
#trp = c("TGG")
#tyr = c("TAC","TAT")
#val = c("GTA","GTT","GTG","GTC")
#----

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
  } else {
    return(c())
  }
}

#----
seqSet = readDNAStringSet("cds/ena-sars.fasta")
seq=seqSet[[11]]
setX = codes.c3[[26]]
setXString = unlist(setX)
setXString = strsplit(setXString, ",")
prefix = setXString[1]
setXString = setdiff(setXString,prefix)
setXString = toString(setXString)
setXString = gsub(", ","",setXString)
setXDNAString = DNAString(setXString)
aa=translate(seq)
#----

codonCount = 0
notPartOfX = 0
replaced = 0

getStatistics = function(){
  return(c(codonCount,notPartOfX,replaced))
}


replaceAminoAcid = function(codonString){
  data("BLOSUM62")
  
}


innerChange = function(codonString){
  
  codonSeq = DNAString(codonString)
  aminoAcid = translate(codonSeq)
  
  codesForThisAA = getCodesForAA(aminoAcid)
  codesForThisAA_size = length(codesForThisAA)
  
  if (!is.null(codesForThisAA)) { #Bei eines Stopp Codon, wird eine leere Liste zurückgegeben
  
    for (j in codesForThisAA_size) {
      testCodon = codesForThisAA[j]
    
      if (grepl(testCodon, setXString)) {
        testCodon = paste("[",testCodon,"]",sep = "")
        replaced<<-replaced+1
       return(testCodon)
      }
    }
    
  }
  
  return(codonString)
}


changeSequence = function (sequence, setX) {
  
  cd=codons(seq) #codons in rna sequence 
  l=length(cd) #amount of codons
  
  codonCount<<-l

  tmpSequence="" #String for sequence after change
  
  for (i in 1:l) {
    
    codon = cd[i]
    codonString = toString(codon)
    
    if (!grepl(codon, setXString)) {
      
      notPartOfX<<-notPartOfX+1

      codonString = innerChange(codonString)
      
    }
    
    #build the new sequence
    tmpSequence = paste(tmpSequence, codonString, sep="")
  }
  return(tmpSequence)
}



#print(toString(seq))
newSequence = changeSequence(seq,setX)
#print(newSequence)









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
  workaround = DNAString("AAA")
  workaround = translate(workaround)
  type1 = typeof(aminoAcid)
  type2 = typeof(workaround)
  if (type1!=type2) {
    aminoAcid = AAString(aminoAcid)
  }
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
# seqSet = readDNAStringSet("cds/ena-sars.fasta")
# seq=seqSet[[11]]
# setX = codes.c3[[26]]
# setXString = unlist(setX)
# setXString = strsplit(setXString, ",")
# prefix = setXString[1]
# setXString = setdiff(setXString,prefix)
# setXString = toString(setXString)
# setXString = gsub(", ","",setXString)
# setXDNAString = DNAString(setXString)
# aa=translate(seq)
#----

codonCount = 0
notPartOfX = 0
replacedCodons = 0
replacedAAs = 0

getStatistics = function(){
  return(c(codonCount,notPartOfX,replacedCodons, replacedAAs))
}


replaceAminoAcid = function(aaString, setXaa, setXString){
  data("BLOSUM62")
  setXaa_length = length(setXaa)
  
  newAA = aaString # sollte sich keine Aminosäure mit score > treshold finden, wird die alte beibehalten
  score = -5 # threshold
  
  for (i in 1:setXaa_length) {
    tmp = BLOSUM62[aaString, setXaa[i]]
    if (tmp >= score) {
      score = tmp
      newAA = setXaa[i]
    }
  }
  if (newAA != aaString) { #wegen treshold notwendig
    replacedAAs<<-replacedAAs+1
  }
  newCodon = innerChange(newAA, setXString)
  return(newCodon)
}


innerChange = function(aminoAcid,setXString,codonString="EMPTY" ){
  
  codesForThisAA = getCodesForAA(aminoAcid)
  codesForThisAA_size = length(codesForThisAA)
  
  if (!is.null(codesForThisAA)) { #Bei eines Stopp Codon, wird eine leere Liste zurückgegeben
  
    for (j in 1:codesForThisAA_size) {
      testCodon = codesForThisAA[j]
    
      if (grepl(testCodon, setXString)) {
        testCodon = paste("[",testCodon,"]",sep = "")
        replacedCodons<<-replacedCodons+1
       return(testCodon)
      }
    }
    
  }
  
  return(codonString)
}


changeSequence = function (sequence, setX) {
  
  cd=codons(seq) #codons in rna sequence 
  l=length(cd) #amount of codons
  
  setXString = unlist(setX)
  setXString = strsplit(setXString, ",")
  prefix = setXString[1]
  setXString = setdiff(setXString,prefix) #Präfix entfernen
  setXString = toString(setXString)
  setXStringPure = gsub(", ","",setXString) #Kommas entfernen
  setXDNAString = DNAString(setXStringPure)
  
  setXaa = translate(setXDNAString)
  setXaa = uniqueLetters(setXaa)
  
  codonCount<<-l

  tmpSequence="" #String for sequence after change
  
  for (i in 1:l) {
    
    codon = cd[i]
    codonString = toString(codon)
    
    if (!grepl(codon, setXString)) {
      
      notPartOfX<<-notPartOfX+1
      
      codonSeq = DNAString(codonString)
      aminoAcid = translate(codonSeq)

      tmp = codonString
      codonString = innerChange(aminoAcid,setXString,codonString)
      
      if (codonString == tmp) {
        codonDNA = DNAString(codonString)
        aaString = toString(translate(codonDNA))
        codonString = replaceAminoAcid(aaString, setXaa, setXString)
      }
      
    }
    
    #build the new sequence
    tmpSequence = paste(tmpSequence, codonString, sep="")
  }
  return(tmpSequence)
}



#print(toString(seq))
#newSequence = changeSequence(seq,setX)
#print(newSequence)









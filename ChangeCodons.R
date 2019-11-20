library(Biostrings)
#source("AminoDecoder.R")

#RNA sequence
seqSet = readDNAStringSet("cds/ena-sars.fasta")
#seqSet = readDNAStringSet("cds/ena-herpes.fasta")
#seqSet = readDNAStringSet("cds/CCDS_nucleotide-human.fasta")
seq=seqSet[[11]]

#Circular Codes
source("ccmotif/R/codes.R")
cCodes = ccmotif.readCodes("ccmotif/C3.txt")
setX = cCodes[23]

setXString = unlist(setX)
setXString = strsplit(setXString, ",")
prefix = setXString[1:2]
setXString = setdiff(setXString,prefix)
setXString = toString(setXString)
setXString = gsub(", ","",setXString)
setXDNAString = DNAString(setXString)

#Amino acids (AA) of interest (AAs coded by non-X codons which can be coded by codons from X))
aaInt=paste("A","F","G","I","L","Q","T","V","Y")
#aaInt = getCodesOfInterest(setX)

#Codons of original sequence
cd=codons(seq)
l=length(cd)

#Amino acid sequence
aa=translate(seq)
print(paste("coding for following amino acids:",aa))



#Counter for statistics
codonCount = l
notPartOfX = 0
changedCodonsCount = 0
mutations = 0
noChangePossibleCount = 0

checkCodon = function(cd){
  
  #String for sequence after change
  tmpSequence=""
  
  for (i in 1:l) {
    #Get codon on position i
    codon = cd[i]
    codonString=toString(codon)
    
    #Trasform Codon into DNAString to translate it
    codonSeq=DNAString(codonString)
    aa1=translate(codonSeq)
    aa1String = toString(aa1)
    
    #Check if codon is not part of set X
    if (!grepl(codonString, setX)) {
      #print(paste(codonString, "ist not part of X. This codon is coding for", aa1String))
      notPartOfX<<-notPartOfX+1
      
      #Check if amino acid (AA) is of interest
      if (aa1String == "*") {
        print("stop-codon")
        noChangePossibleCount<<-noChangePossibleCount+1
      }else if (grepl(aa1String, aaInt)) {
        #print("Gotcha! This amino acid can also be produced by one of the 20 codons in set X")
        changedCodonsCount<<-changedCodonsCount+1
        mutations<<-mutations+1
        
        oldCodon = codonString
        codonString = changeCodon(aa1String,codonString)
        print(paste("Codon changed from",oldCodon, "to", codonString))
        
      }else {
        noChangePossibleCount<<-noChangePossibleCount+1
      }
      
    }
    #build the new sequence
    tmpSequence = paste(tmpSequence, codonString, sep="")
  }
  return(tmpSequence)
}

#depending on the amino acid, replaces the codon with a codon part of X
changeCodon = function(aa,codonString){
  if (aa == "A") {
    codonString="GCC"
  } else if (aa == "F") {
    codonString = "TTC"
  }else if (aa == "G") {
    codonString = "GGC"
  }else if (aa == "I") {
    codonString = "ATC"
  }else if (aa == "L") {
    if(grepl("CT",codonString)){
      codonString = "CTC"
    }else if (grepl("TTA",codonString)) {
      codonString = "CTC"
      #two mutations needed
      mutations<<-mutations+1
    }else if (grepl("TTG",codonString)) {
      codonString = "CTG"
    }
  }else if (aa == "Q") {
    codonString = "CAG"
  }else if (aa == "T") {
    codonString = "ACC"
  }else if (aa == "V") {
    codonString = "GTA"
  }else if (aa == "Y") {
    codonString = "TAC"
  }else{
    print("ERROR")
  }
  return(codonString)
}

newSequenceString = checkCodon(cd)

#new RNA sequence after change
newSequence=DNAString(newSequenceString)
#print(paste("old RNA sequence: ", seq))
#print(paste("new RNA sequence: ", newSequence))

#Amino acid sequence should be the same
newAA=translate(newSequence)
#print(paste("old AA sequence: ", aa))
#print(paste("new AA sequence: ", newAA))
print(paste("The AA sequence is the same:", aa==newAA))

print(paste("The sequence consists of", codonCount, "codons"))
print(paste(notPartOfX,"(",floor((notPartOfX/codonCount)*100), "% ) of these codons are NOT part of set X"))
print(paste(changedCodonsCount, "(", floor((changedCodonsCount/codonCount)*100), "% ) have been changed through", mutations, "mutations without altering the amino acid sequence"))
print(paste(noChangePossibleCount,"(",floor((noChangePossibleCount/codonCount)*100), "% ) could not be changed without altering the amino acid sequence"))

#TODO: genauen Grad der Ver?nderung in % angeben (wie viele Basen insgesamt ge?ndert wurden)

print("END")




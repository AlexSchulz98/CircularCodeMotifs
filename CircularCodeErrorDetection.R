#20 Codons which are part of the circular code
circularCode=c("AAC", "AAT", "ACC", "ATC", "ATT", "CAG", "CTC", "CTG", "GAA", "GAC",  "GAG", "GAT", "GCC", "GGC", "GGT", "GTA", "GTC", "GTT", "TAC", "TTC")

#Counter for statistics
combinations = 0
afterOne = 0
afterTwo = 0
afterThree = 0
afterFour = 0


for (t in 1:20000) {
  ran = floor(runif(5, min=1, max=20))
  
  codon1= circularCode[ran[1]]
  codon2= circularCode[ran[2]]
  codon3= circularCode[ran[3]]
  codon4= circularCode[ran[4]]
  codon5= circularCode[ran[5]]
  
  seq = paste(codon1, codon2, codon3, codon4, codon5, sep="")
  
  #frame shift by 1
  seqfs1 = sub(".","",seq)
  
  #frame shift by 2
  seqfs2 = sub(".","",seqfs1)
  
  #get new codons for each frame shift
  newseq1 = sapply(seq(from=1, to=nchar(seqfs1), by=3), function(i) substr(seqfs1, i, i+2))
  newseq2 = sapply(seq(from=1, to=nchar(seqfs2), by=3), function(i) substr(seqfs2, i, i+2))
  
  testCodons = function(fs, newseq){
    combinations<<-combinations+1
    for(x in 1:20){
      if (grepl(newseq[1], circularCode[x], fixed=TRUE)) {
        #print(paste("Erstes Codon bei Frame Shift", fs, "ergibt noch Sinn"))
        #print(newseq)
        afterOne<<-afterOne+1
        for(y in 1:20){
          if(grepl(newseq[2], circularCode[y], fixed=TRUE)) {
            print(paste("1.+ 2. Codon bei Frame Shift", fs, "ergibt noch Sinn. Selten"))
            print(newseq)
            afterTwo<<-afterTwo+1
            for(z in 1:20){
              if(grepl(newseq[3], circularCode[z], fixed=TRUE)) {
                print(paste("1.+ 2.+ 3. Codon bei Frame Shift", fs, "ergibt noch Sinn. Außergewöhnlich"))
                print(newseq)
                afterThree<<-afterThree+1
                for(zz in 1:20){
                  if(grepl(newseq[4], circularCode[zz], fixed=TRUE)) {
                    print(paste("1.+ 2.+ 3.+ 4. Codon bei Frame Shift", fs, "ergibt noch Sinn. Das hätte nicht passieren dürfen"))
                    print(newseq)
                    afterFour<<-afterFour+1
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  #frameshift = 1
  testCodons(1, newseq1)
  
  #frameshift = 2
  testCodons(2, newseq2)
}

print(paste("Kombinationen insgesamt:",combinations))
print(paste("Bei",afterOne, "Kombinationen war das erste Codon noch Teil von X"))
print(paste("Bei",afterTwo, "Kombinationen war das zweite Codon noch Teil von X"))
print(paste("Bei",afterThree, "Kombinationen war das dritte Codon noch Teil von X"))
print(paste("Bei",afterFour, "Kombinationen war das vierte Codon noch Teil von X"))











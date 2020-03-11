library(Biostrings)
library(ccmotif)
library(seqinr)
source("../Scripts/Mutations.R")
source("../Scripts/Sequences.R")


seqSet_input = readDNAStringSet("cds_ex/ena-sars.fasta") # change here
seqName = "SarsExample" # change here
circularCodes = c(23) # change here
frame = 0 # change here
subAminos = FALSE # change here
threshold = 0 # change here

seqSet = seqSet_input

if (frame > 0) {
  #change reading frame
  for (j in 1:length(seqSet_input)) {
    seqSet[[j]] = changeReadingFrame(frame, seqSet_input[[j]])
  }
}

# mutate sequence for each circular code
for (h in 1:length(circularCodes)) {
  ccNumber = circularCodes[h]
  setX = codes.c3[[ccNumber]]
  
  # mutate every sequence in this set
  for (z in 1:length(seqSet)) {
    newSeq = mutateSequence(seqSet[[z]], setX, subAminos, threshold)
    
    # write result into fasta file
    write.fasta(
      sequences = newSeq,
      names = paste("modified", names(seqSet[z])),
      file.out = paste(seqName,
                       "_",
                       ccNumber,
                       "_",
                       frame,
                       ".fasta",
                       sep = ""),
      open = "a",
      nbchar = 60,
      as.string = FALSE
    )
    
  }
  
}
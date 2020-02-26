library("Biostrings")
library(seqinr)
source("Scripts/Sequences.R")

#' Set path
#' Set file
#' Set frame
#' Output: Fasta file
#' deletes the first [frame] bases from each cds in the set


path ="../BA Circular Code/cds/" #change here
file ="orf_genomic_yeast.fasta"  #change here
frame = 1 # change here


seqSet=readDNAStringSet(paste(path,file,sep=""))

for (i in 1:length(seqSet)) { 
  seqSet[[i]] = changeReadingFrame(frame, seqSet[[i]])
  
  write.fasta(
    sequences = seqSet[i],
    names = paste("frameshifted", names(seqSet[i])),
    file.out = paste(
      "cds/",
      "frame_",
      frame,
      "_",
      file,
      sep = ""
    ),
    open = "a",
    nbchar = 60,
    as.string = FALSE
  )
}

















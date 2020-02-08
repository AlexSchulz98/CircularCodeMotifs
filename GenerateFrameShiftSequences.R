library("ccmotif") # Version 0.6.6
library("Biostrings")
library(seqinr)
source("RawDataExtraction.R")
source("CodeManipulation.R")

#' Input is a fasta file (seqSet) and the wanted frame
#' deletes the first [frame] bases from each sequences in the set
#' output is a fasta file
#' care for length of new file, currently capped to 1000


path ="../BA Circular Code/cds/"
file ="orf_genomic_yeast.fasta"  #change here
frame = 2 # change here

print(file)
print(frame)


seqSet=readDNAStringSet(paste(path,file,sep=""))

for (i in 1:pmin(length(seqSet), 1000)) { #change length if wanted
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

















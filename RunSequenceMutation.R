library(Biostrings)
library(ccmotif)
library(seqinr)
source("CodeManipulation.R")
source("Main.R")

#' input is a fasta file [seqSet], the sequence name [seqName], a list of circular codes [circularCodes] and the frame [frame]
#' at max the first [setLength] sequences in a set are used
#' output is a fasta file (of max. [setLength] sequences)

start_time_global <- Sys.time()

seqSet_input = readDNAStringSet("cds/ena-sars.fasta") # change here
seqName = "SarsVirus" # change here          
circularCodes = c(188,51,49) # change here
frame = 1 # change here
setLength = 1000 #change here

seqSet = seqSet_input
  
  if (frame > 0) {
    #change reading frame
    for (j in 1:pmin(length(seqSet_input), setLength)) {
      seqSet[[j]] = changeReadingFrame(frame, seqSet_input[[j]])
    }
  }
  
  for (h in 1:length(circularCodes)) {
    
    start_time_code <- Sys.time()
    
    ccNumber = circularCodes[h] #circular codes for each organism specified in Parameter.R
    setX = codes.c3[[ccNumber]]
    
    
    for (z in 1:pmin(length(seqSet), setLength)) {
      # first 1000 sequences are examined
      
      start_time_change <- Sys.time()
      
      newSeq = main(seqSet[[z]], setX)

      write.fasta(
        sequences = newSeq,
        names = paste("modified", names(seqSet[z])),
        file.out = paste(
          "output_sequences_new/",
          seqName,
          "_",
          ccNumber, 
          "_",
          frame,
          ".fasta",
          sep = ""
        ),
        open = "a",
        nbchar = 60,
        as.string = FALSE
      )
      
      end_time_change <- Sys.time()
      print(paste(ccNumber,"- Done with set", z, "/", pmin(length(seqSet), setLength), "in"))
      print(end_time_change - start_time_change)
    }
    
    end_time_code <- Sys.time()
    print(paste("fasta file generated:", seqName,"code",ccNumber,"frame",frame))
    print(end_time_code - start_time_code)
    
  }
  
end_time_global <- Sys.time()
print("Finished")
print(end_time_global - start_time_global)







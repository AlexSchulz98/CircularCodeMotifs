source("Parameter.R")

start_time_global <- Sys.time()

seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
seqName = "SarsVirus" #for naming generated files
circularCodes = c(84,82,206,23)

for (i in 1) {
  start_time_frame <- Sys.time()
  
  if (i > 0) {
    #change reading frame
    for (j in 1:pmin(length(seqSet), 1000)) {
      seqSet[[j]] = changeReadingFrame(i, seqSet[[j]])
    }
  }
  
  for (h in 1:length(circularCodes)) {
    
    start_time_code <- Sys.time()
    
    ccNumber = circularCodes[h] #circular codes for each organism specified in Parameter.R
    setX = codes.c3[[ccNumber]]
    
    
    for (z in 1:pmin(length(seqSet), 1000)) {
      # first 1000 sequences are examined
      
      start_time_change <- Sys.time()
      
      newSeq = main(seqSet[[z]], setX)

      write.fasta(
        sequences = newSeq,
        names = paste("modified", names(seqSet[z])),
        file.out = paste(
          "output_sequences/output_",
          seqName,
          "_",
          setX[[1]],
          "_frame_",
          i,
          ".fasta",
          sep = ""
        ),
        open = "a",
        nbchar = 60,
        as.string = FALSE
      )
      
      end_time_change <- Sys.time()
      print(paste("Done with set", z, "/", pmin(length(seqSet), 1000), "in"))
      print(end_time_change - start_time_change)
    }
    
    end_time_code <- Sys.time()
    print(paste("fasta file generated:", seqName,"code",setX[[1]],"frame",i))
    print(end_time_code - start_time_code)
    
  }
  
  end_time_frame <- Sys.time()
  print(paste("Done with frame", i))
  print(end_time_frame - start_time_frame)
  
}

# saveRDS(object = outputMatrix_codons, file = paste("output_matrixes_codons/codons_outputMatrix_",seqName,".RDS", sep = ""))
# write.xlsx(outputMatrix_codons, file ="output_matrixes_codons/codons_outputMatrix.xlsx", sheetName = seqName,col.names = TRUE, row.names = TRUE, append = TRUE)
# 
# saveRDS(object = outputMatrix_aminos, file = paste("output_matrixes_aminos/aminos_outputMatrix_",seqName,".RDS", sep = ""))
# write.xlsx(outputMatrix_aminos, file ="output_matrixes_aminos/aminos_outputMatrix.xlsx", sheetName = seqName,col.names = TRUE, row.names = TRUE, append = TRUE)

end_time_global <- Sys.time()
print("Finished")
print(end_time_global - start_time_global)







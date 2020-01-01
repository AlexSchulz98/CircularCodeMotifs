source("Parameter.R")

start_time_global <- Sys.time()

for (i in 1:length(codes.c3)) {
  
  start_time_code <- Sys.time()
  
  setX = codes.c3[[i]]
  
  for (z in 1:length(seqSet)) {
    
    start_time_change <- Sys.time()
    
    newSeq = main(seqSet[[z]], setX)
    
    write.fasta(
      sequences = newSeq,
      names = paste("modified", names(seqSet[z])),
      file.out = paste("output_sequences/output_",seqName,"_", setX[[1]], ".fasta", sep = ""),
      open = "a",
      nbchar = 60,
      as.string = FALSE
    )
    
    # outputMatrix_codons = codonCount(codon_table, length(codes.c3), seqSet[[z]], newSeq, i)
    # codon_table = outputMatrix_codons
    # 
    # outputMatrix_amino = aminoCount(amino_table, length(codes.c3), seqSet[[z]], newSeq, i)
    # amino_table = outputMatrix_amino
    
    end_time_change <- Sys.time()
    print(paste("Done with set", z, "/", length(seqSet),"in"))
    print(end_time_change - start_time_change)
  }
  
  end_time_code <- Sys.time()
  print(paste("Done with circular code", i, "/", length(codes.c3)))
  print(end_time_code - start_time_code)
  
}

# saveRDS(object = outputMatrix_codons, file = paste("output_matrixes_codons/codons_outputMatrix_",seqName,".RDS", sep = ""))
# write.xlsx(outputMatrix_codons, file ="output_matrixes_codons/codons_outputMatrix.xlsx", sheetName = seqName,col.names = TRUE, row.names = TRUE, append = TRUE)
# 
# saveRDS(object = outputMatrix_aminos, file = paste("output_matrixes_aminos/aminos_outputMatrix_",seqName,".RDS", sep = ""))
# write.xlsx(outputMatrix_aminos, file ="output_matrixes_aminos/aminos_outputMatrix.xlsx", sheetName = seqName,col.names = TRUE, row.names = TRUE, append = TRUE)

end_time_global <- Sys.time()
print("Finished")
print(end_time_global - start_time_global)







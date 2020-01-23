
# -----
# for (i in 1:3) { #length(codes.c3)
#   
#   start_time_run <- Sys.time()
#   
#   fileName = paste("output_sequences/output_SarsVirus_X", i, ".fasta", sep ="")
#   newSeq = readDNAStringSet(fileName)
#   
#   for (z in 1:length(seqSet)) {
#     start_time_tables <- Sys.time()
#     
#     outputMatrix_codons = codonCount(codon_table, length(codes.c3), seqSet[[z]], newSeq[[z]], i)
#     codon_table = outputMatrix_codons
#     
#     
#     outputMatrix_aminos = aminoCount(amino_table, length(codes.c3), seqSet[[z]], newSeq[[z]], i)
#     amino_table = outputMatrix_aminos
#     
#     end_time_tables <- Sys.time()
#     
#     print(paste("Done with set",z,"/",length(seqSet),"(X-Code:",i,")"))
#     print(end_time_tables - start_time_tables)
#     
#   }
#   
#   end_time_run <- Sys.time()
#   print(paste("X-Code",i,"finished"))
#   print(end_time_run - start_time_run)
# }
# ----

# ###### Frameshift beachten! #######
# dnaf = readDNAStringSet(file)
# 
# for (j in 1:pmin(length(dnaf), 1000)) {
#   dnaf[[j]] = changeReadingFrame(1, dnaf[[j]])
# }

start_time_files <- Sys.time()

ar = generateEmptyTable(64,64,CODONS)

dnafmod = readDNAStringSet(file23_2) #change here
seqName = "Herpes_X23_frame_2" #change here

for (i in 1:length(dnafmod)) {
  outputMatrix_codons = codonCount(ar,dnaf[[i]], dnafmod[[i]])
  ar = outputMatrix_codons
}

saveRDS(
  object = outputMatrix_codons,
  file = paste(
    "output_matrixes_codons/codons_outputMatrix_",
    seqName,
    ".RDS",
    sep = ""
  )
)

end_time_files <- Sys.time()
print(end_time_files - start_time_files)

# write.xlsx(
#   rds23_1,
#   file = paste("output_matrixes_codons/codons_outputMatrix_",
#                seqName,
#                ".xlsx",
#                sep=""),
#   sheetName = "Sheet1",
#   col.names = TRUE,
#   row.names = TRUE,
#   append = TRUE
# )


#----
# # RDS codons
# saveRDS(
#   object = outputMatrix_codons,
#   file = paste(
#     "output_matrixes_codons/codons_outputMatrix_",
#     seqName,
#     ".RDS",
#     sep = ""
#   )
# )
# 
# #RDS amino acids
# saveRDS(
#   object = outputMatrix_aminos,
#   file = paste(
#     "output_matrixes_aminos/aminos_outputMatrix_",
#     seqName,
#     ".RDS",
#     sep = ""
#   )
# )

# write.xlsx(
#   outputMatrix_codons,
#   file = "output_matrixes_codons/codons_outputMatrix.xlsx",
#   sheetName = seqName,
#   col.names = TRUE,
#   row.names = TRUE,
#   append = TRUE
# )
# 
# write.xlsx(
#   outputMatrix_aminos,
#   file = "output_matrixes_aminos/aminos_outputMatrix.xlsx",
#   sheetName = seqName,
#   col.names = TRUE,
#   row.names = TRUE,
#   append = TRUE
# )




#Nützliche CCMotif Funktionen:
# ccmotif.scan.fasta --> Analysis of motif lengths in a list of fasta sequences (spuckt 3 datein aus)





source("Parameter.R")
source("RawDataExtraction.R")

path = "D:/Projekte/BA Circular Code/output_sequences/EscherichiaColi/" #change here
fastafile = list.files(path, pattern = "*.fasta") 

dnaf = readDNAStringSet("cds/Escherichia_coli.HUSEC2011CHR1.cds.all.fasta") #change here

dnaf1 = dnaf
for (j in 1:pmin(length(dnaf), 1000)) {
  dnaf1[[j]] = changeReadingFrame(1, dnaf[[j]])
}
dnaf2 = dnaf
for (j in 1:pmin(length(dnaf), 1000)) {
  dnaf2[[j]] = changeReadingFrame(2, dnaf[[j]])
}


for (h in 1:length(fastafile)) {
  
  print(fastafile[h])
  tmp = unlist(strsplit(fastafile[h],"_"))
  seqName = tmp[1]
  code = tmp[2]
  frame = unlist(strsplit(tmp[3],".fasta"))
  
  ar = generateEmptyTable(64, 64, CODONS)
  
  dnafmod = readDNAStringSet(paste(path,fastafile[h],sep=""))
  
  if (as.numeric(frame) == 0) {
    for (i in 1:length(dnafmod)) {
      outputMatrix_codons = codonCount(ar, dnaf[[i]], dnafmod[[i]])
      ar = outputMatrix_codons
      print(paste("Done with sequence", i, "/",length(dnafmod)))
    }
  } else if (as.numeric(frame) == 1) {
    for (i in 1:length(dnafmod)) {
      outputMatrix_codons = codonCount(ar, dnaf1[[i]], dnafmod[[i]])
      ar = outputMatrix_codons
      print(paste("Done with sequence", i, "/",length(dnafmod)))
    }
  } else if (as.numeric(frame) == 2) {
    for (i in 1:length(dnafmod)) {
      outputMatrix_codons = codonCount(ar, dnaf2[[i]], dnafmod[[i]])
      ar = outputMatrix_codons
      print(paste("Done with sequence", i, "/",length(dnafmod)))
    }
  } else {
    print("Error")
  }
  
  saveRDS(
    object = outputMatrix_codons,
    file = paste(
      "output_matrixes_codons/",
      seqName,
      "_",
      code,
      "_",
      frame,
      ".RDS",
      sep = ""
    )
  )
  
}


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


#Nützliche CCMotif Funktionen:
# ccmotif.scan.fasta --> Analysis of motif lengths in a list of fasta sequences (spuckt 3 datein aus)





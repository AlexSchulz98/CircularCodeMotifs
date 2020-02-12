library(Biostrings)
source("RawDataExtraction.R")
source("CodeManipulation.R")
source("Parameter.R")

#' Input: directory only with mutated fasta files of the original sequence
#' reads in all files and generates a 64x64 for each

path = "../BA Circular Code/output_sequences_new/Rein/" # change mutated sequence set here
fastafile = list.files(path, pattern = "*.fasta") 

dnaf = readDNAStringSet("cds/ena-ch-reinhardtii.fasta") #change original sequence set here

dnaf = dnaf[1:1009] # change here for deleting IUPAC Codes
dnaf = deleteIUPACSequences(dnaf)

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


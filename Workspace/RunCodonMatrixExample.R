library(Biostrings)
source("../Scripts/DataAnalysis.R")
source("../Scripts/Sequences.R")

dnaf = readDNAStringSet("cds_ex/ena-sars.fasta") # original sequence

path = "../Workspace/" # change here
fastafile = list.files(path, pattern = "*.fasta")

for (h in 1:length(fastafile)) {
  dnafmod = readDNAStringSet(paste(path, fastafile[h], sep = "")) # modified sequence
  
  tmp = unlist(strsplit(fastafile[h], "_"))
  seqName = tmp[1]
  code = tmp[2]
  frame = unlist(strsplit(tmp[3], ".fasta"))
  
  ar = generateEmptyTable(64, 64, CODONS) # new matrix
  
  # fill matrix
  for (i in 1:length(dnafmod)) {
    outputMatrix_codons = codonCount(ar, dnaf[[i]], dnafmod[[i]])
    ar = outputMatrix_codons
    print(paste("Done with sequence", i, "/", length(dnafmod)))
  }
  
  # generate RDS file
  saveRDS(object = outputMatrix_codons,
          file = paste(seqName,
                       "_",
                       code,
                       "_",
                       frame,
                       ".RDS",
                       sep = ""))
}

library(Biostrings)
library(ccmotif)
library(xlsx)
source("CodeManipulation.R")
source("Main.R")
source("RawDataExtraction.R")

# 64 codons
CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
# 20 amino acids + 1 *  for stop codons
AMINOACIDS = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
#setX = codes.c3[[26]] #Circular Code

codon_table = generateEmptyTable(64,64,length(codes.c3), CODONS) #empty table for codons
amino_table = generateEmptyTable(21,21,length(codes.c3), AMINOACIDS) #empty table for amino acisd

#---------------------------------------------

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
      file.out = paste("output_sequences/output_SarsVirus_", setX[[1]], ".fasta", sep = ""),
      open = "a",
      nbchar = 60,
      as.string = FALSE
    )
    
    outputMatrix_codons = codonCount(codon_table, length(codes.c3), seqSet[[z]], newSeq, i)
    codon_table = outputMatrix_codons
    
    outputMatrix_amino = aminoCount(amino_table, length(codes.c3), seqSet[[z]], newSeq, i)
    amino_table = outputMatrix_amino
    
    end_time_change <- Sys.time()
    print(paste("Done with set", z, "/", length(seqSet),"in",end_time_change - start_time_change))
  }
  
  end_time_code <- Sys.time()
  print(paste("Done with circular code", i, "/", length(codes.c3),"in",end_time_code - start_time_code))
  
}

saveRDS(object = outputMatrix_codons, file = "output_matrixes_codons/codons_outputMatrix_SarsVirus.RDS")
write.xlsx(outputMatrix_codons, file ="output_matrixes_codons/codons_outputMatrix.xlsx", sheetName = "Sars Virus",col.names = TRUE, row.names = TRUE, append = TRUE)

saveRDS(object = outputMatrix_aminos, file = "output_matrixes_aminos/aminos_outputMatrix_SarsVirus.RDS")
write.xlsx(outputMatrix_aminos, file ="output_matrixes_aminos/aminos_outputMatrix.xlsx", sheetName = "Sars Virus",col.names = TRUE, row.names = TRUE, append = TRUE)

end_time_global <- Sys.time()
print(paste("Finished in",end_time_global - start_time_global))








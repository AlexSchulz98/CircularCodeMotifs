library(Biostrings)
library(ccmotif)
library(xlsx)
source("CodeManipulation.R")
source("Main.R")
source("RawDataExtraction.R")

#----

# 64 codons
CODONS = c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT","CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC" ,"GTA" ,"GTG", "GCT" ,"GCC" ,"GCA" ,"GCG" ,"GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
# 20 amino acids + 1 *  for stop codons
AMINOACIDS = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

#----


seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
seqName = "SarsVirus" #for naming generated files

codon_table = generateEmptyTable(64,64,length(codes.c3), CODONS) #empty table for codons
amino_table = generateEmptyTable(21,21,length(codes.c3), AMINOACIDS) #empty table for amino acisd
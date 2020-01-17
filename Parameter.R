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
circularCodes = c(3,23,41,161)

# seqSet = readDNAStringSet("cds/CCDS_nucleotide-human.fasta") # DNA (RNA) sequence set
# seqName = "Human" #for naming generated files
# circularCodes = c(23,25,117,173)

# seqSet = readDNAStringSet("cds/celegans.fasta") # DNA (RNA) sequence set
# seqName = "Celegans" #for naming generated files
# circularCodes = c(23,24,122,171)

# seqSet = readDNAStringSet("cds/ena-ch-reinhardtii.fasta") # DNA (RNA) sequence set
# seqName = "Reinhardtii" #for naming generated files
# circularCodes = c(10,23,166,173)

# seqSet = readDNAStringSet("cds/ena-herpes.fasta") # DNA (RNA) sequence set
# seqName = "Herpes" #for naming generated files
# circularCodes = c(23,172,173)

# seqSet = readDNAStringSet("cds/orf_genomic_yeast.fasta") # DNA (RNA) sequence set
# seqName = "Yeast" #for naming generated files
# circularCodes = c(20,23,122,171)

# seqSet = readDNAStringSet("cds/Escherichia_coli.HUSEC2011CHR1.cds.all.fasta")
# seqName = "EscherichiaColi"
# circularCodes = c(20,22,23)


codon_table = generateEmptyTable(64,64,length(codes.c3), CODONS) #empty table for codons
amino_table = generateEmptyTable(21,21,length(codes.c3), AMINOACIDS) #empty table for amino acisd
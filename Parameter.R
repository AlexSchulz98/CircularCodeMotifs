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


#161,41,3
# seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
# seqName = "SarsVirus" #for naming generated files

#25,173,117
seqSet = readDNAStringSet("cds/CCDS_nucleotide-human.fasta") # DNA (RNA) sequence set
seqName = "Human" #for naming generated files

#122,24,171
# seqSet = readDNAStringSet("cds/celegans.fasta") # DNA (RNA) sequence set
# seqName = "Celegans" #for naming generated files

#166,173,10
# seqSet = readDNAStringSet("cds/ena-ch-reinhardtii.fasta") # DNA (RNA) sequence set
# seqName = "Reinhardtii" #for naming generated files

# 173,172,23
# seqSet = readDNAStringSet("cds/ena-herpes.fasta") # DNA (RNA) sequence set
# seqName = "Herpes" #for naming generated files

# 122,171,20
# seqSet = readDNAStringSet("cds/orf_genomic_yeast.fasta") # DNA (RNA) sequence set
# seqName = "Yeast" #for naming generated files

# 23,20,22
# seqSet = readDNAStringSet("cds/Escherichia_coli.HUSEC2011CHR1.cds.all.fasta")
# seqName = "EscherichiaColi"

# ???
# seqSet = readDNAStringSet("cds/CCDS_nucleotide.mouse.current.fasta") 
# seqName = "Mouse" 


codon_table = generateEmptyTable(64,64,length(codes.c3), CODONS) #empty table for codons
amino_table = generateEmptyTable(21,21,length(codes.c3), AMINOACIDS) #empty table for amino acisd
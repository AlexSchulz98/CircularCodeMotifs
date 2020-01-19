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

# seqSet = readDNAStringSet("cds/ena-sars.fasta") # DNA (RNA) sequence set
# seqName = "SarsVirus" #for naming generated files
# circularCodes0 = c(3,23,41,161)
# circularCodes1 = c(84,82,206,23)
# circularCodes2 = c(59,55,193,23)

# seqSet = readDNAStringSet("cds/CCDS_nucleotide-human.fasta") # DNA (RNA) sequence set
# seqName = "Human" #for naming generated files
# circularCodes0 = c(23,25,117,173)
# circularCodes1 = c(52,131,55,23)
# circularCodes2 = c(59,88,69,23)

# WARNING: Invalide Basen
# seqSet = readDNAStringSet("cds/celegans.fasta") # DNA (RNA) sequence set
# seqName = "Celegans" #for naming generated files
# circularCodes0 = c(23,24,122,171)
# circularCodes1 = c(73,95,76,23)
# circularCodes2 = c(193,192,88,23)

# WARNING: Invalide Basen
# seqSet = readDNAStringSet("cds/ena-ch-reinhardtii.fasta") # DNA (RNA) sequence set
# seqName = "Reinhardtii" #for naming generated files
# circularCodes0 = c(10,23,166,173)
# circularCodes1 = c(191,187,194,23)
# circularCodes2 = c(156,150,100,23)

# seqSet = readDNAStringSet("cds/ena-herpes.fasta") # DNA (RNA) sequence set
# seqName = "Herpes" #for naming generated files
# circularCodes0 = c(23,172,173)
# circularCodes1 = c(191,42,126,23)
# circularCodes2 = c(210,89,88,23)

# seqSet = readDNAStringSet("cds/orf_genomic_yeast.fasta") # DNA (RNA) sequence set
# seqName = "Yeast" #for naming generated files
# circularCodes0 = c(20,23,122,171)
# circularCodes1 = c(141,70,143,23)
# circularCodes2 = c(100,65,150,23)

# seqSet = readDNAStringSet("cds/Escherichia_coli.HUSEC2011CHR1.cds.all.fasta")
# seqName = "EscherichiaColi"
# circularCodes0 = c(20,22,23)
# circularCodes1 = c(86,90,87,23)
# circularCodes2 = c(189,193,61,23)


codon_table = generateEmptyTable(64,64, CODONS) #empty table for codons
amino_table = generateEmptyTable(21,21, AMINOACIDS) #empty table for amino acisd
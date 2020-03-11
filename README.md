# On mutations in coding sequences with circular code motifs

Abstract:
Circular Codes offer a possible solution to explain the high error tolerance of the genetic code. By using circular codes a frameshift can be detected within two to threecodons. Therefore a continuous sequence of circular code codons, a so called motif,is essential. The ancient genetic code is believed to consist only of codons from circular codes. To test this theory, a selection of coding sequences was mutated based on 216 self-complementary C 3 codes. The aim was to restore the ancient sequence, considering conditions for a biologically realistic substitution. Sequences with a code usage of up to 84.32 % could be generated. In two sequences the average length of the motifs could be increased by a total of 25.32% compared to the expected value. However, the best results were achieved outside of the original reading frame. Therefore a further tracking of shifted sequences is suggested.

R version: 3.6.2 (2019-12-12)
Libraries: Biostrings, seqinr, ccmotif (M. Gumbel, version 0.6.6), ggplot2, xlsx, knitr

../cds: coding sequences from different organisms

../codes: self-complementary-c3-codes

../Scripts: core functions
Mutations.R: methods used to mutate sequences
DataAnalysis.R: methods for metrics and 64x64 codon matrix
Sequences.R: manipulation of sequences (frameshift, UIPAC ambiguity codes, etc)

../Workspace: example
1. Start RunMutationExample.R
--> generates mutated sequence (fasta file)
2. Start RunCodonMatrixExample.R
--> generates codon matric (rds file)
3. Start Example_Analysis
--> see results as dataframe and graph


How to use main program: 

- set correct paths in scripts
- modify marked variables in script
- generate mutated sequence (fasta file) via RunSequenceMutation.R
- generate codon matrix (rds file) via RunCodonSequence.R
- run and see results in Overall_Analysis.Rmd
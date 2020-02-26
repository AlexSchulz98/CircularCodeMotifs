library("ccmotif") # Version 0.6.6

#' Set path
#' Set file
#' Output: RDS
#' calculates Code Usage for all 216 codes and saves data to file

path = "../BA Circular Code/" # Change here
file = "Sars_23_0.fasta" # Change here


ff = read.fasta(paste(path, file, sep = ""),
                as.string = TRUE,
                forceDNAtolower = FALSE)


codons = codon.splitlist(ff)
cu = codon.usage(codons)

filename = unlist(strsplit(file, ".fasta"))

df = data.frame(File = c(),
                Code = c(),
                CodeUsage = c())

for (h in 1:length(codes.c3)) {
  code = h
  p = codes.usage(cu, codes.c3[[code]]) #Code Usage
  
  
  #Fill data frame
  tmp = data.frame(File = file,
                  Code = code,
                  CodeUsage = round(p, 3))
  
  df = rbind(df, tmp)
  
}

#generate RDS file
subFrame = df[df$File == file,]
subFrame = subFrame[order(-subFrame$CodeUsage),] #sort

saveRDS(
  object = subFrame,
  file = paste("CodeUsage/",
               "CU",
               "_",
               filename,
               ".RDS",
               sep = "")
)
  



















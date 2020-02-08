library("ccmotif") # Version 0.6.6
source("RawDataExtraction.R")
source("CodeManipulation.R")

#' input is a fasta file
#' care for right frame
#' output is an RDS file 
#' containing all stop-codon-free c3-codes + the corresponding code usage in % 
#' order is descending

path = "../BA Circular Code/cds/" # (Change here)
file = "frame_2_orf_genomic_yeast.fasta" # Change here

df = data.frame(File=c(),
                Code=c(),
                CodeUsage=c())

  
  ff = read.fasta(paste(path,file,sep=""), as.string = TRUE, forceDNAtolower = FALSE)
  
  print(file)
  
  codons = codon.splitlist(ff)
  cu = codon.usage(codons)
  
  for (h in 1:length(codes.c3)) {
    
    stop = codes.containsStop(codes.c3[[h]]$codons)
    
    if (stop == TRUE) {
      
      code = h
      print(code)
      
      p = codes.usage(cu, codes.c3[[code]]) #Code Usage
      
      print(paste(code, p))
      
      #Fills data frame
      tmp = data.frame(File=file,
                       Code=code,
                       CodeUsage=round(p,3))
      df = rbind(df,tmp)
      
    }
  }

  #generates RDS file
  print(file)
  
  subFrame = df[df$File==file,]
  subFrame = subFrame[order(-subFrame$CodeUsage),]
  
  saveRDS(
    object = subFrame,
    file = paste(
      "CodeUsage/",
      "CodeUsage",
      "_",
      file,
      "_",
      "NoStopCodon",
      ".RDS",
      sep = ""
    )
  )
  



















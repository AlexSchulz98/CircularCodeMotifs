# Adhoc analysis
library(ggplot2)
#library(ccmotif)
source("analysis.R")

E = function(p) 1 / (1 - p)
V = function(p) p / (1 - p)^2

ccmotif.adhoc = function(filename, codes = codes.c3[23],
                         frame = 0,
                         minCu = .6, maxCu = .8,
                         tmpDir = "adhoc/") {

  s = list(id = filename, fid = filename, idxRange = c(1, 1))

  doScan = function(filename, codes, frame = 0) {
    filepath = paste(tmpDir, filename, ".fasta", sep = "")
    print(filepath)
    ff = read.fasta(filepath, as.string = TRUE,
                    forceDNAtolower = FALSE)
    label = paste(s$fid, "_1_1", sep = "")
    ccmotif.scan.fasta(ff, codes = codes, frame = frame,
                       label = label, tmpDir = tmpDir)
  }

  # TODO Redundant to noteobook.
  plotCuMlRnd = function(s, frame = 0, motif = TRUE,
                         minCu = .2, maxCu = .9) {
    label = if (motif) {
      paste(s, " in frame ", frame, sep = "")
    } else {
      paste(s, " (non motif) in frame ", frame, sep = "")
    }
    col = if (motif) 2 else 3

    cu = codes.readCodeUsage(s, tmpDir,
                             frame, rnd = FALSE)
    ml = ccmotif.readMotifLengths(s, tmpDir,
                                  frame, rnd = FALSE, motif = motif,
                                  newZipVersion = TRUE)

    mml = ccmotif.motiflengthAggr(ml, mean, "mean")
    M = cbind(mml, cu = cu[, col]) # code usage and mean motif lengths

    cuv = seq(minCu, maxCu, length.out = 20)
    mlv = E(cuv)
    vlv = mlv + sqrt(V(cuv))
    df = data.frame(M)
    ggplot() +
      geom_point(data = df,
                 aes(x = cu, y = mean, color = "C3 codes")) +
      geom_line(data = data.frame(x = cuv, y = mlv),
                aes(x = x, y = y, color = "expected")) +
      # geom_line(data = data.frame(x = cuv, y = vlv),
      #          aes(x = x, y = y, color = "stdev")) +
      scale_color_manual(values = c("C3 codes" = "blue",
                                    "random codes" = "red",
                                    "expected" = "black",
                                    "stdev" = "gray")) +
      labs(x = "code usage", y = "average motif length",
           title = label) +
      theme_bw()
  }

  doScan(filename, codes, frame = frame)
  p = plotCuMlRnd(s, frame = frame, minCu = minCu, maxCu = maxCu)
  plot(p)
}


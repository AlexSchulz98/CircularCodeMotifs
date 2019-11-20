# Routines for analysis.

# Workaround. Set this global variable:
# newZipVersion = TRUE

ccmotif.motiflengthAggr = function(ml, op, label) {
  mml = apply(ml, 2, op, na.rm = TRUE) # aggregate column wise
  dmml = cbind.data.frame(as.character(colnames(ml)), mml)
  colnames(dmml) = c("code", label)
  dmml
}


#' Chi square test for distribution. (Bosch, page 379)
#'
#' Might be replaced by a pure R routine - if useful?
#'
#' @param ml List of motif lengths
#' @param p Probability (relative frequency) of code usage
#' @param a alpha
#'
#' @return FALSE if hypothesis is rejected, TRUE if not.
ccmotif.chi2motiflengths = function(ml, p, a) {
  n = length(na.omit(ml)) # Number of motif lengths
  cl = ccmotif.classes(ml, p)
  N = cl$sample; n0 = cl$geom
  m = 1:length(N)
  h = sapply(m, function(j) (N[j] - n0[j])^2 / n0[j])
  chi2 = sum(sapply(m, function(j) (N[j] - n0[j])^2 / n0[j]))
  # chi2 <= qchisq(1 - a, df = n - 1)
  pchisq(chi2, df = n - 1, lower.tail = FALSE)
}

# global variable
mlCache = list()

#' Read and cache motif lengths.
#'
#' @param s
#' @param results_path_prefix
#' @param frame
#' @param rnd
#' @param motif
#'
#' @return
ccmotif.cacheMotifLengths = function(s, results_path_prefix,
                                     frame = 0, rnd = FALSE, motif = TRUE) {
  # Calculate unique file name.
  # Warning: redundant to ccmotif.readMotifLengths()
  frameLabel = paste("_f", frame, sep = "")
  rndLabel = if (rnd) "_rnd" else ""
  motifLabel = if (motif) "" else "non_"
  rangeIdx = paste("_", s$idxRange[1], "_", s$idxRange[2],
                   sep = "")
  mlFile = paste("motif_lengths_", motifLabel, s$fid, rndLabel, rangeIdx,
                 frameLabel, ".csv", sep = "")

  l = length(mlCache) # Number of cached lists.
  if (is.null(mlCache[[mlFile]])) { # not yet cached?
    ml = ccmotif.readMotifLengths(s, results_path_prefix, frame, 
                                  rnd, motif, 
                                  newZipVersion = newZipVersion)
    if (l <= 30) { # Limit cache size to avoid memory problems.
      mlCache[[mlFile]] <<- ml # global var!
      if (mlCacheInfo) {
        print(paste(mlFile, " cached. Size: ", length(mlCache), sep = ""))
      }
    }
  } else {
    ml = mlCache[[mlFile]]
    if (mlCacheInfo) {
      print(paste(mlFile, " from cache.", sep = ""))
    }
  }
  ml # return motif lengths
}


# equivalence classes

ccmotif.ec.avgCodeUsage = function(species_settings) {
  # frame with dummy values:
  df = data.frame(c(list(""), 1:27), stringsAsFactors = FALSE)
  colnames(df) = c("Species", 1:27)
  for (s in species_settings) {
    cu = codes.readCodeUsage(s, results_path_prefix)
    m = unlist(lapply(1:27, function(i) mean(cu[codes.c3.inclass(i), 2])))
    tmp = data.frame(c(list(s$id), m), stringsAsFactors = FALSE)
    colnames(tmp) = c("Species", 1:27)
    df = rbind(df, tmp)
  }
  df = data.frame(df[2:dim(df)[1], ])
  colnames(df) = c("Species", 1:27)
  df
}


ccmotif.ec.avgMotifLength = function(species_settings, ec) {
  clidx = codes.c3.inclass(ec)
  M = matrix(clidx, nrow = 1, byrow = TRUE)

  for (s in species_settings) {
    ml = ccmotif.cacheMotifLengths(s, results_path_prefix)
    mml = ccmotif.motiflengthAggr(ml, mean, "mean")
    m = unlist(mml[clidx, 2])
    M = rbind(M, m)
  }
  M = M[2:(length(species_settings) + 1), ]
  colnames(M) = clidx
  rownames(M) = unlist(lapply(species_settings, function(s) s$id))

  A = matrix(apply(M, 2, mean), nrow = 1, byrow = TRUE)
  colnames(A) = clidx
  B = A[, order(A[1, ], decreasing = TRUE)]
  B
}

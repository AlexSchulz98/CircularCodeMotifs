
foo = function() {
  amcode = c("AAC", "AAT", "ACC", "ATC", "ATT", "CAG", "CTC", "CTG", "GAA", "GAC", "GAG", "GAT", "GCC", "GGC", "GGT", "GTA", "GTC", "GTT", "TAC", "TTC")
  i = 1
  for (code in codes.c3) {
    if (identical(amcode, code)) {
      print(code)
      print(i)
    }
    i = i + 1
  }
}

b2v = function(b) {
  if (b == "A") return(0)
  if (b == "T") return(1)
  if (b == "C") return(2)
  if (b == "G") return(3)
}

c2v = function(c) {
  s = strsplit(c, split = "")[[1]]
  16 * b2v(s[1]) + 4 * b2v(s[2]) + b2v(s[3])
}
trim <- function(x, ctrim = c.trim, ntrim = n.trim) {
  substring(x, ctrim, nchar(x)-ntrim)
}

Jaccard = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

getTCR <- function(tmp) {
  chains <- list(one = c("TCR1","cdr3_aa1"), two = c("TCR2","cdr3_aa2"))
  TCR <- NULL
  for (i in 1:2) {
    sub <- as.data.frame(tmp[,c("barcode", chains[[i]])])
    colnames(sub) <- c("barcode", "TCR", "cdr3_aa")
    
    sub$v <- str_split(sub$TCR, "[.]", simplify = T)[,1]
    sub$j <- str_split(sub$TCR, "[.]", simplify = T)[,2]
    colnames(sub)[3] <- "Var1"
    TCR[[i]] <- sub
    sub <- NULL
  }
  return(TCR)
}
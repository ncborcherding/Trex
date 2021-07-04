#method is the distance calculation
#group is how to recombine the list before the calculation
#group should be placed in list, such as group = list(set1 = c("names", "list" "elements"), set2 = etc)
#Threshold for nearest neighbor identification
#matrix is the type of matrix to return for the calculation 
#TRA - for normalized TCRA, TRB for normalized TCRB, MED - mean edit distance for both loci, 
#nn - nearest neighbor, or jaccard for similarity index matrix
#c.trim - number of AA to trim from the front of the cdr3 sequence
#n.trim - number of AA to trim from the end of the cdr3 sequence
#barcodes - corresponding cell names from the single cell object used to filter

library(stringdist)
distanceMatrix <- function(TCR, 
                           edit.method = "lv",
                           threshold = threshold, 
                           c.trim = c.trim,
                           n.trim = n.trim) {
knn.return <- list()
length.check <- min(nchar(TCR[[i]]$cdr3_aa))
  for (i in seq_along(TCR)) {
    if (c.trim == 0 & n.trim == 0){
      TR <- TCR[[i]]$cdr3_aa
    } else {
      TR<- trim(TCR[[i]]$cdr3_aa, c.trim, n.trim)
    }
  barcodes <- TCR[[i]]$barcode
  length <- nchar(as.character(TR))
  length.check <- min(length)
  if (length.check < n.trim + c.trim + 3) {
    print("Current trim strategy leaves less than 3 AA residues calucations, please consider 
        prefilter short cdr3 aa sequences
        or reset the trimming parameters")
    if (length.check < n.trim + c.trim) {
      stop("Unable to perform edit distance calculations, a cdr3 AA sequence is shorted than
           the trimming parameters")
    }
  }
  TR <- as.matrix(stringdistmatrix(TR, method = edit.method))
  #This converts the distance matrices calculated above to a normalized value based on the length of the cdr3 sequence.
  medianlength <- median(na.omit(length))
  out_matrix <- matrix(ncol = ncol(TR), nrow=nrow(TR))
    for (k in seq_len(ncol(TR))) {
      for (l in seq_len(nrow(TR))) {
        if (is.na(length[l]) | is.na(length[k])) {
          out_matrix[k,l] <- NA
          out_matrix[l,k] <- NA
        } else {
          if (length[k] - length[l] >= round(medianlength/1.5)) {
            out_matrix[k,l] <- 1 - (TR[k,l]/(max(length[k], length[l])))
            out_matrix[l,k] <- 1 - (TR[l,k]/(max(length[k], length[l])))
          } else {
            out_matrix[k,l] <- 1 - (TR[k,l]/((length[k]+ length[l])/2))
            out_matrix[l,k] <- 1 - (TR[l,k]/((length[k]+ length[l])/2))
          }
      }
      }
    }
    out_matrix  <- TCR[[i]]$barcode
    out_matrix <- TCR[[i]]$barcode
  #Converting relative edit distance to adjacency matrix
  out_matrix[out_matrix >= threshold] <- 1
  
  knn.matrix <- matrix(ncol=ncol(out_matrix), nrow=nrow(out_matrix), 0)
  for (m in 1:nrow(knn.matrix)){
    # find closes neighbors - not technically nearest neighbor
    matches <- which(out_matrix[m,] > threshold) #all neighbors with > 0.8 edit similarity
    knn.matrix[m,matches] <- 1
    knn.matrix[matches,m] <- 1
  } 
  rownames(knn.matrix) <- TCR[[i]]$barcode
  colnames(knn.matrix) <- TCR[[i]]$barcode
  knn.return[[i]] <- knn.matrix
  }
return(knn.return)
}


#Chains for MAIT cells
#' @importFrom dplyr bind_rows
scoreMAIT <- function(TCR, species = NULL) {
  membership <- bind_rows(TCR)
  comp <- list(mouse = list(v = "TRAV1", j = "TRAJ33", length = 12), 
               human = list(v = "TRAV1-2", j = c("TRAJ33", "TRAJ20", "TRAJ12"), length = 12))
  score <- data.frame("barcode" = unique(membership[,"barcode"]), score = 0)
  cells <- unique(score$barcode)
  for (i in seq_len(length(cells))) {
    v <- membership[membership[,1] == cells[i],]$v
    j <- membership[membership[,1] == cells[i],]$j
    length <- nchar(membership[membership[,1] == cells[i],]$Var1)
    if(comp[[species]]$v  %in% v && comp[[species]]$j %in% j && comp[[species]]$length %in% length) {
      score$score[i] <- 1
    } else {
      next()
    }
  }
  return(score)
}
#Chains for INKT cells
#' @importFrom dplyr bind_rows
scoreINKT <- function(TCR, species = NULL) {
  membership <- bind_rows(TCR)
  comp <- list(mouse = list(v = "TRAV11", j = "TRAJ18", length = 15), 
               human = list(v = "TRAV10", j = c("TRAJ18", "TRBV25"), length = c(14,15,16)))
  score <- data.frame("barcode" = unique(membership[,"barcode"]), score = 0)
  cells <- unique(score$barcode)
  for (i in seq_len(length(cells))) {
    v <- membership[membership[,1] == cells[i],]$v
    j <- membership[membership[,1] == cells[i],]$j
    length <- nchar(membership[membership[,1] == cells[i],]$Var1)
    if(comp[[species]]$v  %in% v && comp[[species]]$j %in% j && comp[[species]]$length %in% length) {
      score$score[i] <- 1
    } else {
      next()
    }
  }
  return(score)
}

#results in 1- Relative distance between V and J regions of TCRA
alphaDistance <- function(getTCR) {
  references <- read.delim("./data/imgt_tra_locus_order.txt")
  references$IMGT.gene.name <-  stringr::str_remove_all(references$IMGT.gene.name, "/")
  
  membership <- getTCR[[1]]
  score <- data.frame("barcode" = membership[,"barcode"], score = 0)
  cells <- unique(score$barcode)
  for (i in seq_len(length(cells))) {
    v.ref <- references[,1][which(references[,1] == membership[membership$barcode == cells[i], ]$v)]
    j.ref <- references[,1][which(references[,1] == membership[membership$barcode == cells[i], ]$j)]
    if (length(v.ref) > 0 & length(v.ref) > 0) {
      score$score[i] <- 1-(references[references[,1] == j.ref,]$IMGT.gene.order - 
                             references[references[,1] == v.ref,]$IMGT.gene.order) / length(references[,2])
    } else {
      score$score[i] <- NA
    }
  }
  return(score)
}

#results in 1- Relative distance between V and J regions of TCRB 
#To Do: Think about influence of direction on the calculation
betaDistance <- function(getTCR) {
  references <- read.delim("./data/imgt_trb_locus_order.txt", header = FALSE)
  colnames(references) <- c("IMGT.gene.name", "IMGT.gene.order", "direction")
  references$IMGT.gene.name <-  stringr::str_remove_all(references$IMGT.gene.name, "/")
  
  membership <- getTCR[[2]]
  score <- data.frame("barcode" = membership[,"barcode"], score = 0)
  cells <- unique(score$barcode)
  for (i in seq_len(length(cells))) {
    v.ref <- references[,1][which(references[,1] == membership[membership$barcode == cells[i], ]$v)]
    j.ref <- references[,1][which(references[,1] == membership[membership$barcode == cells[i], ]$j)]
    if (length(v.ref) > 0 & length(v.ref) > 0) {
      score$score[i] <- 1-(references[references[,1] == j.ref,]$IMGT.gene.order - 
                             references[references[,1] == v.ref,]$IMGT.gene.order) / length(references[,2])
    } else {
      score$score[i] <- NA
    }
  }
  return(score)
}


aaProperty <- function(TCR, 
                       c.trim = 0,
                       n.trim = 0) { 
  aa.score <- list()
  reference <- read.delim("./data/aa_props.tsv")
  for (i in seq_along(TCR)) {
    membership <- TCR[[i]]
    score <- as.data.frame(matrix(ncol = 29, nrow = length(unique(membership[,"barcode"]))))
    colnames(score) <- c("barcodes", colnames(reference)[2:29])
    score$barcodes <- unique(membership[,"barcode"])
    cells <- unique(score$barcode)
    for (j in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[j],]$Var1
      if (c.trim != 0 | n.trim != 0){
        tmp.CDR <- trim(tmp.CDR, ctrim = c.trim, ntrim = n.trim)
      }
      refer <- unlist(strsplit(tmp.CDR, ""))
      int <- reference[match(refer, reference$aa),]
      score[j,2:29] <- colSums(int[,2:29])/length(refer)
      #score[j,2:29] <- colMedians(as.matrix(int[,2:29]))
    } 
    aa.score[[i]] <- score
  }
  return(aa.score)
}

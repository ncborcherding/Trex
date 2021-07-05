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
distanceMatrix <- function(TCR.list, 
                           edit.method = "lv",
                           edit.threshold = threshold, 
                           ctrim = c.trim,
                           ntrim = n.trim) {
  knn.return <- list()
  
    for (i in seq_along(TCR.list)) {
      if (ctrim == 0 & ntrim == 0){
        TR <- TCR.list[[i]]$cdr3_aa
      } else {
        TR<- trim(TCR.list[[i]]$cdr3_aa, ctrim, ntrim)
      }
    barcodes <- TCR.list[[i]]$barcode
    length <- nchar(na.omit(as.character(TR)))
    length.check <- min(length)
    if (length.check < ntrim + ctrim + 3 && ntrim + ctrim > 0) {
      print("Current trim strategy leaves less than 3 AA residues calculations, please consider
          prefilter short cdr3 aa sequences or changing the trimming parameters")
      if (length.check < ntrim + ctrim) {
        stop("Unable to perform edit distance calculations, atleast one cdr3 AA sequence is
        shorter than the trimming parameters")
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
    
    knn.matrix <- matrix(ncol = ncol(TR), nrow=nrow(TR))
    for (m in 1:nrow(knn.matrix)){
      # find closes neighbors - not technically nearest neighbor
      matches <- which(out_matrix[m,] > edit.threshold) #all neighbors with > 0.8 edit similarity
      knn.matrix[m,matches] <- 1
      knn.matrix[matches,m] <- 1
    } 
    rownames(knn.matrix) <- TCR.list[[i]]$barcode
    colnames(knn.matrix) <- TCR.list[[i]]$barcode
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
                       c.trim = c.trim,
                       n.trim = n.trim, 
                       AA.properties = AA.properties) { 
  aa.score <- list()
  
  reference <- read.delim("./data/aa_props.tsv")
  col.ref <- grep(tolower(paste(AA.properties, collapse = "|")), colnames(reference))
  other.ref <- grep("af|kf", colnames(reference)[-1], invert = TRUE)
  if ("other" %in% AA.properties | AA.properties == "all") {
    column.ref <- sort(c(col.ref, other.ref))
  } else {
    column.ref <- sort(col.ref)
  }
  for (i in seq_along(TCR)) {
    membership <- TCR[[i]]
    score <- as.data.frame(matrix(ncol = length(column.ref)+1, nrow = length(unique(membership[,"barcode"]))))
    colnames(score) <- c("barcodes", colnames(reference)[column.ref])
    score$barcodes <- unique(membership[,"barcode"])
    cells <- unique(score$barcode)
    for (j in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[j],]$Var1
      if (c.trim != 0 | n.trim != 0){
        tmp.CDR <- trim(tmp.CDR, ctrim = c.trim, ntrim = n.trim)
      }
      refer <- unlist(strsplit(tmp.CDR, ""))
      int <- reference[match(refer, reference$aa),]
      score[j,col.ref] <- colSums(int[,col.ref])/length(refer)
    } 
    dist <- as.matrix(Dist(score[,seq_len(ncol(score))[-1]], method = "pearson"))
    max <- max(dist)
    dist <- (max-dist)/max
    knn.matrix <- matrix(ncol=ncol(dist), nrow=nrow(dist), 0)
    rownames(knn.matrix) <- names
    colnames(knn.matrix) <- names
    for (m in 1:nrow(dist)) {
      matches <- which(dist[m,] > threshold)
      knn.matrix[m,matches] <- 1
      knn.matrix[matches,m] <- 1
    }
    aa.score[[i]] <- knn.matrix
  }
  return(aa.score)
}

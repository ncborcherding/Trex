#Calculating Edit Distance Matrix
#' @importFrom stringdist stringdistmatrix
distanceMatrix <- function(TCR, 
                           edit.method = "lv",
                           nearest.method = nearest.method,
                           near.neighbor = near.neighbor,
                           threshold = threshold, 
                           n.trim = n.trim,
                           c.trim = c.trim)
 {
    return <- list()
    for (i in seq_along(TCR)) {
      if (c.trim != 0 || n.trim != 0 ){
        TR<- trim(TCR[[i]]$cdr3_aa, c.trim, n.trim)
      } else {
        TR <- TCR[[i]]$cdr3_aa
      }
      barcodes <- TCR[[i]]$barcode
      length <- nchar(na.omit(as.character(TR)))
      length.check <- min(length)
      if (length.check < (n.trim + c.trim + 3) && c.trim != 0 && n.trim != 0) { 
        warning(strwrap(prefix = " ", initial = "", "Current trim strategy leaves less 
        than 3 AA residues calculations, please consider
            prefiltering short cdr3 AA sequences or changing the trimming parameters"))}
      if (length.check < (n.trim + c.trim)) {
          stop(strwrap(prefix = " ", initial = "", "Unable to perform edit distance 
          calculations, at least one cdr3 AA sequence is
          shorter than the trimming parameters"))
      }
    TR <- as.matrix(stringdistmatrix(TR, method = edit.method))
    #This converts the distance matrices calculated above to a normalized 
    #value based on the length of the cdr3 sequence.
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
    return[[i]] <- out_matrix
    }
  names(return) <- paste0(names(TCR), ".edit")
  return(return)
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

#Calculating Distance of AA in CDR3 using mean
#' @importFrom amap Dist
aaProperty <- function(TCR, 
                       c.trim = c.trim,
                       n.trim = n.trim, 
                       nearest.method = nearest.method,
                       near.neighbor = near.neighbor,
                       threshold = threshold,
                       AA.method = AA.method,
                       AA.properties = AA.properties) { 
  aa.score <- list()
  load("./data/Trex.Data.rda") ### Need to add reference data
  reference <- Trex.Data[[1]]
  col.ref <- grep(tolower(paste(AA.properties, collapse = "|")), colnames(reference))
  if (AA.properties == " both") {
    column.ref <- unique(sort(c(AF.col, KF.col)))
  } else {
    column.ref <- unique(sort(col.ref))
  }
  chain <- names(TCR)
  for (i in seq_along(TCR)) {
    membership <- TCR[[i]]
    names <- membership$barcode
    if (AA.method != "auto") {
      score <- as.data.frame(matrix(ncol = length(column.ref)+1, nrow = length(unique(membership[,"barcode"]))))
      colnames(score) <- c("barcodes", colnames(reference)[column.ref])
      score$barcodes <- unique(membership[,"barcode"])
      array.reshape <- NULL
    else {
      aa.model <- aa.model.loader(chain[[i]], AA.properties)
      range <- aa.range.loader(chain[[i]], AA.properties, Trex.Data) 
      local.min <- range[[1]]
      local.max <- range[[2]]
    }
    cells <- unique(membership[,"barcode"])
    for (j in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[j],]$cdr3_aa
      if (AA.method != "auto") {
        if (c.trim != 0 | n.trim != 0){
          tmp.CDR <- trim(tmp.CDR, c.trim = c.trim, n.trim = n.trim)
        }
        refer <- unlist(strsplit(tmp.CDR, ""))
        int <- reference[match(refer, reference$aa),]
        score[j,column.ref] <- colSums(int[,column.ref])/length(refer)
      } else {
        refer <- unlist(strsplit(tmp.CDR, ""))
        refer <- c(refer, rep(NA, 50 - length(refer)))
        int <- reference[match(refer, reference$aa),c(1,col.ref)]
        array.reshape.tmp <- array_reshape(as.matrix(int[,-1]), length(col.ref)*50)
        array.reshape <- rbind(array.reshape, array.reshape.tmp)
        next()
      } 
    }
    if (aa.method == "auto") {
      #Here is where the autoencoder embeds and returns a 30-vector value for each cdr3
      score <- auto.embedder(array.reshape, aa.model, local.max, local.min)
      score <- data.frame(unique(membership[,"barcode"]), score)
    }
    dist <- as.matrix(Dist(score[,seq_len(ncol(score))[-1]], method = "pearson"))
    max <- max(dist, na.rm = TRUE)
    dist <- (max-dist)/max
    aa.score[[i]] <- dist
  }
  return(aa.score)
}

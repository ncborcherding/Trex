trim <- function(x, ctrim = c.trim, ntrim = n.trim) {
  substring(x, ctrim, nchar(x)-ntrim)
}

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

Jaccard = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

gene.to.knn <- function(tmpscore) {
  names <- tmpscore$barcode
  knn.matrix <- matrix(ncol=length(names), nrow=length(names), 0)
  rownames(knn.matrix) <- names
  colnames(knn.matrix) <- names
  match <- which(tmpscore > 0)
  knn.matrix[match,match] <- 1
  return(knn.matrix)
}

'%!in%' <- Negate("%in%")
add.to.network <- function(network, new.knn, name) {
  length <- length(network)
  names.list <- names(network)
  for (i in seq_along(list(new.knn))) {
    network[length + i] <- new.knn
    names(network) <- c(names.list, name)
  }
  return(network)
}

#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment coldata 'coldata<-'
add.meta.data <- function(sc, meta) {
if (inherits(x=sc, what ="Seurat")) { 
  col.name <- names(meta) %||% colnames(meta)
  sc[[col.name]] <- meta
} else {
  rownames <- rownames(colData(sc))
  colData(sc) <- cbind(colData(sc), meta[rownames,])[, union(colnames(colData(sc)),  colnames(meta))]
  rownames(colData(sc)) <- rownames  
}
  return(sc)
}
  
#' Function to pull and organize TCR depending on the chain selected
#' @importFrom stringr str_split
getTCR <- function(sc, chains) {
  meta <- grabMeta(sc)
  tmp <- data.frame(barcode = rownames(meta), 
                    str_split(meta[,"CTaa"], "_", simplify = TRUE), 
                    str_split(meta[,"CTgene"], "_", simplify = TRUE))
  if (chains %in% c("TRA", "TRD")) {
    pos <- list(c(2,4))
  } else if (chains %in% c("TRB", "TRG")) {
    pos <- list(c(3,5))
  } else {
    pos <- list(one = c(2,4), two = c(3,5))
    ch.1 <- grep("TRB|TRA",sc[[]]$CTgene[1])
    if (ch.1 == 1) {
      chains <- c("TRA", "TRB")
    } else {
      chains <- c("TRD", "TRG")
    }
  }
  TCR <- NULL
  for (i in seq_along(pos)) {
    sub <- as.data.frame(tmp[,c(1,pos[[i]])])
    colnames(sub) <- c("barcode", "cdr3_aa", "genes")
    sub$v <- str_split(sub$genes, "[.]", simplify = T)[,1]
    sub$j <- str_split(sub$genes, "[.]", simplify = T)[,2]
    sub[sub == ""] <- NA
    TCR[[i]] <- sub
    sub <- NULL
  }
  names(TCR) <- chains
  return(TCR)
}

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData 
grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  }
  else if (inherits(x=sc, what ="SummarizedExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}

#This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}


library(muxViz)
multiplex.network <- function(multi.network) {
  layers <- length(multi.network)
  layerCouplingStrength <- 1
  networkOfLayersType <- "categorical"
  isDirected <- FALSE
  nodeTensor <- list() 
  g.list <- list() 
  for (l in 1:layers) {
    #Generate the layers
    g.list[[l]] <- graph_from_adjacency_matrix(master.list[[l]], mode = "undirected")
    #Get the list of adjacency matrices which build the multiplex
    nodeTensor[[l]] <- get.adjacency(g.list[[l]])
  }
  layerTensor <- BuildLayersTensor(Layers=Llyers, OmegaParameter=layerCouplingStrength,
                                   MultisliceType=networkOfLayersType)
  M <- BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(nodeTensor, layerTensor, layers, Nodes)
  N <- GetAggregateNetworkFromSupraAdjacencyMatrix(M, layers, Nodes)
  N <- get.adjacency(N)
}







#T
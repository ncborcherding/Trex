
#Allows for the cutting of the cdr3 sequenc
trim <- function(x, ctrim = c.trim, ntrim = n.trim) {
  substring(x, ctrim, nchar(x)-ntrim)
}

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

#Convert binary to adjacency matrix
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
#Append the multiple network object
add.to.network <- function(network, new.knn, name) {
  length <- length(network)
  names.list <- names(network)
  if(is(new.knn)[1] != "list") { list(new.knn) }
  for (i in seq_along(new.knn)) {
    network[[length + i]] <- new.knn[[i]]
    
  }
  names(network) <- c(names.list, name)
  return(network)
}

# Add to meta data some of the metrics calculated
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
  if (length(chains) == 1 & chains != "both") {
    if (chains %in% c("TRA", "TRD")) {
      pos <- list(c(2,4))
    } else if (chains %in% c("TRB", "TRG")) {
      pos <- list(c(3,5))
    }
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
  else if (inherits(x=sc, what ="SingleCellExperiment")){
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

#Shhhhhh
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#This multiplexes the network and gets simplified eigen values
#' @importFrom muxViz BuildLayersTensor BuildSupraAdjacencyMatrixFromEdgeColoredMatrices GetAggregateNetworkFromSupraAdjacencyMatrix
#' @importFrom igraph simplify spectrum
multiplex.network <- function(multi.network, n.dim, sc) {
  Nodes <- nrow(multi.network[[1]])
  layers <- length(multi.network)
  layerCouplingStrength <- 1
  networkOfLayersType <- "categorical"
  isDirected <- FALSE
  nodeTensor <- list() 
  g.list <- list() 
  for (l in seq_len(layers)) {
    #Generate the layers
    g.list[[l]] <- graph_from_adjacency_matrix(multi.network[[l]], mode = "undirected")
    #Get the list of adjacency matrices which build the multiplex
    nodeTensor[[l]] <- get.adjacency(g.list[[l]])
  }
  layerTensor <- BuildLayersTensor(Layers=layers, OmegaParameter=layerCouplingStrength,
                                   MultisliceType=networkOfLayersType)
  M <- BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(nodeTensor, layerTensor, layers, Nodes)
  N <- GetAggregateNetworkFromSupraAdjacencyMatrix(M, layers, Nodes)
  N <- simplify(N)
  eigen <- spectrum(N, 
                    which = list(howmany = n.dim), 
                    algorithm = "arpack")
  eigen <- eigen$vectors
  rownames(eigen) <- rownames(sc[[]])
  colnames(eigen) <- paste0("Trex_", seq_len(ncol(eigen)))
  return(eigen)
}

#Define adjacency matrix by either threshold or nearest neighbor
#' @importFrom FNN knn.index
get.knn <- function(TCR, i, nearest.method, near.neighbor, edit.threshold) {
  knn.matrix <- matrix(ncol = nrow(TCR[[i]]), nrow= nrow(TCR[[i]]))
  if (nearest.method == "threshold") {
    for (m in seq_len(nrow(knn.matrix))){
      # find closes neighbors - not technically nearest neighbor
      matches <- which(out_matrix[m,] > edit.threshold) #all neighbors > threshold similarity
      knn.matrix[m,matches] <- 1
      knn.matrix[matches,m] <- 1
    }
    }else if (nearest.method == "nn") {
      matches <- knn.index(out_matrix,k = near.neighbor)
      for (m in seq_len(nrow(matches))) {
        neigh.check <- which(out_matrix[m,] == 1) 
        if (length(neigh.check) > near.neighbor) {
          matches <- sample(neigh.check, near.neighbor)
        }
        knn.matrix[m,matches] <- 1
        knn.matrix[matches,m] <- 1
      }
    } 
    rownames(knn.matrix) <- TCR[[i]]$barcode
    colnames(knn.matrix) <- TCR[[i]]$barcode
    return(knn.matrix)
}

#Add the eigen values to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
adding.DR <- function(sc, maTrex, reduction.name) {
  if (inherits(sc, "Seurat")) {
    DR <- CreateDimReducObject(
      embeddings = maTrex,
      loadings = maTrex,
      projected = maTrex,
      stdev = rep(0, ncol(maTrex)),
      key = reduction.name,
      jackstraw = NULL,
      misc = list())
    sc[[reduction.name]] <- DR
  } else if (inherits(sc, "SingleCellExperiment")) {
    reducedDim(sc, reduction.name) <- maTrex
  }
  
  
}

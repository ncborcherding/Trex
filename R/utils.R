
# Allows for the cutting of the cdr3 sequence
trim <- function(x, n.trim = n.trim, c.trim = c.trim) {
  substring(x, n.trim, nchar(x)-c.trim)
}

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

# Convert binary to adjacency matrix
gene.to.knn <- function(tmpscore) {
  names <- tmpscore$barcode
  knn.matrix <- matrix(ncol=length(names), nrow=length(names), 0)
  rownames(knn.matrix) <- names
  colnames(knn.matrix) <- names
  match <- which(tmpscore > 0)
  knn.matrix[match,match] <- 1
  return(knn.matrix)
}
#Invert %in%
'%!in%' <- Negate("%in%")

# Append the multiple network object
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
#' @importFrom SingleCellExperiment colData
add.meta.data <- function(sc, meta, header) {
  barcodes <- meta$barcode
  meta <- as.data.frame(meta[,2])
  colnames(meta) <- header
  rownames(meta) <- barcodes
  
if (inherits(x=sc, what ="Seurat")) { 
  col.name <- names(meta) %||% colnames(meta)
  sc[[col.name]] <- meta
} else {
  rownames <- rownames(colData(sc))
  colData(sc) <- cbind(colData(sc), 
          meta[rownames,])[, union(colnames(colData(sc)),  colnames(meta))]
  rownames(colData(sc)) <- rownames  
}
  return(sc)
}
  
#Function to pull and organize TCR depending on the chain selected
#' @importFrom stringr str_split
getTCR <- function(sc, chains) {
  meta <- grabMeta(sc)
  tmp <- data.frame(barcode = rownames(meta), 
                    str_split(meta[,"CTaa"], "_", simplify = TRUE), 
                    str_split(meta[,"CTgene"], "_", simplify = TRUE))
  if (length(chains) == 1 && chains != "both") {
    if (chains %in% c("TRA", "TRD")) {
      pos <- list(c(2,4))
    } else if (chains %in% c("TRB", "TRG")) {
      pos <- list(c(3,5))
    }
  } else {
    pos <- list(one = c(2,4), two = c(3,5))
    ch.1 <- grep("TRB|TRA",sc[[]]$CTgene[1])
    chains <- c("TRA", "TRB")
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
#' @importFrom igraph simplify spectrum graph_from_adjacency_matrix get.adjacency
multiplex.network <- function(multi.network, n.dim, barcodes) {
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
                    which = list(howmany = n.dim+1), 
                    algorithm = "arpack")
  eigen <- eigen$vectors[,seq_len(n.dim)]
  rownames(eigen) <- barcodes
  colnames(eigen) <- paste0("Trex_", seq_len(ncol(eigen)))
  return(eigen)
}


#Retunrs appropriate model for autoencoder
#' @importFrom tensorflow tf
#' @importFrom keras load_model_hdf5
aa.model.loader <- function(chain, AA.properties) {
  quiet(tensorflow::tf$compat$v1$disable_eager_execution())
    select  <- system.file("extdata", paste0(chain, "_", 
                               AA.properties, "_Encoder.h5"), 
                          package = "Trex")
  model <- quiet(load_model_hdf5(select, compile = FALSE))
  return(model)
}

#Selects columns to normalize input data basedon the inputs to the model
aa.range.loader <- function(chain, AA.properties, Trex.Data) {
  range <- Trex.Data[["model.ranges"]][[chain]]
  min <- range[["min"]]
  max <- range[["max"]]
  ref <- seq(1, 750, 15)
  if (AA.properties == "AF") {
    ref2 <- sort(c(ref, ref+1, ref+2, ref+3, ref+4))
    min <- min[ref2]
    max <- max[ref2]
  } else if (AA.properties == "AF") {
    ref2 <- sort(c(ref+5, ref+6, ref+7, ref+8, ref+9, ref+10, ref+11, ref+12, ref+13, ref+14))
    min <- min[ref2]
    max <- max[ref2]
  }
  range <- list(min = min, max = max)
  return(range)
}




#Define adjacency matrix by either threshold or nearest neighbor
#' @importFrom FNN knn.index
get.knn <- function(barcodes, out_matrix, nearest.method, near.neighbor, threshold) {
  if (nearest.method == "threshold") {
    knn.matrix <- matrix(ncol = nrow(out_matrix), nrow= nrow(out_matrix))
    for (m in seq_len(nrow(knn.matrix))){
      # find closes neighbors - not technically nearest neighbor
      matches <- which(out_matrix[m,] > threshold) #all neighbors > threshold similarity
      knn.matrix[m,matches] <- 1
      knn.matrix[matches,m] <- 1
    }
  } else if (nearest.method == "nn") {
      knn.matrix <- matrix(ncol = nrow(out_matrix), nrow= nrow(out_matrix))
      knn <- knn.index(out_matrix,k = near.neighbor)
      for (m in seq_len(nrow(knn.matrix))) {
        neigh.check <- which(out_matrix[m,] == 1) 
        if (length(neigh.check) > near.neighbor) {
          matches <- sample(neigh.check, near.neighbor)
        } else {
          matches <- knn[m,]
        }
        knn.matrix[m,matches] <- 1
        knn.matrix[matches,m] <- 1
      }
    } 
    rownames(knn.matrix) <- barcodes
    colnames(knn.matrix) <- barcodes
    return(knn.matrix)
}

#Add the eigen vectors to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
adding.DR <- function(sc, reduction, reduction.name) {
  if (inherits(sc, "Seurat")) {
    DR <- suppressWarnings(CreateDimReducObject(
      embeddings = reduction,
      loadings = reduction,
      projected = reduction,
      stdev = rep(0, ncol(reduction)),
      key = reduction.name,
      jackstraw = NULL,
      misc = list()))
    sc[[reduction.name]] <- DR
  } else if (inherits(sc, "SingleCellExperiment")) {
    reducedDim(sc, reduction.name) <- reduction
  }
  return(sc)
  
}

AF.col <- c(2,3,4,5,6)
KF.col <- c(7,8,9,10,11,12,13,14,15,16)

#Generats the 30 vector based on autoencoder model 
#First normalizes the value by the min and max of the autoencoder training data
auto.embedder <- function(array.reshape, aa.model, local.max, local.min) {
  for(i in seq_len(ncol(array.reshape))) {
    (array.reshape[,i] - local.min[i])/(local.max[i] - local.min[i])
  }
  array.reshape[is.na(array.reshape)] <- 0
  score <- stats::predict(aa.model, array.reshape)
  return(score)
}
  


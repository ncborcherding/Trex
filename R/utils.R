
# Allows for the cutting of the cdr3 sequence
trim <- function(x, n.trim = n.trim, c.trim = c.trim) {
  substring(x, n.trim, nchar(x)-c.trim)
}

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

# Convert binary to adjacency matrix
#' @importFrom igraph graph_from_adjacency_matrix
gene.to.knn <- function(tmpscore) {
  names <- tmpscore$barcode
  match <- which(tmpscore$score > 0)
  knn.norm = list(data.frame("from" = match[1],
                          "to" = match))
  return(knn.norm)
}
#Invert %in%
'%!in%' <- Negate("%in%")


# Add to meta data some of the metrics calculated
#' @importFrom rlang %||%
#' @importFrom SingleCellExperiment colData
add.meta.data <- function(sc, meta, header) {
  barcodes <- rownames(meta)
  
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
    sub$v <- str_split(sub$genes, "[.]", simplify = TRUE)[,1]
    sub$j <- str_split(sub$genes, "[.]", simplify = TRUE)[,2]
    sub[sub == ""] <- NA
    TCR[[i]] <- sub
    sub <- NULL
  }
  names(TCR) <- chains
  return(TCR)
}

#This is to grab the meta data from a Seurat or SCE object
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

#This multiplexes the network and gets simplified eigenvectors
#' @importFrom muxViz BuildLayersTensor BuildSupraAdjacencyMatrixFromExtendedEdgelist GetAggregateNetworkFromSupraAdjacencyMatrix
#' @importFrom igraph simplify spectrum graph_from_edgelist get.edgelist cluster_louvain
multiplex.network <- function(multi.network, n.dim, barcodes) {
  Nodes <- length(barcodes)
  layers <- length(multi.network)
  networkOfLayersType <- "categorical"
  nodeTensor <- NULL 
  g.list <- list() 
  for (l in seq_len(layers)) {
    #Get the list of adjacency matrices which build the multiplex
    g.list[[l]] <- graph_from_edgelist(as.matrix(multi.network[[l]]), directed = FALSE)
    tmp <-  get.edgelist(g.list[[l]])
    tmp <- data.frame("node.from" = tmp[,1], "layer.from" = l, "node.to" = tmp[,2], "layer.to" = l,  "weight" = 1)
    nodeTensor <- rbind.data.frame(nodeTensor, tmp)
  }
  
  layerTensor <- BuildLayersTensor(Layers=layers, OmegaParameter=1,
                                   MultisliceType=networkOfLayersType)
  M <- suppressMessages(BuildSupraAdjacencyMatrixFromExtendedEdgelist(nodeTensor, layers, Nodes, isDirected = FALSE))
  N <- GetAggregateNetworkFromSupraAdjacencyMatrix(M, layers, Nodes)
  N <- simplify(N)
  eigen <- spectrum(N, 
                    which = list(howmany = n.dim), 
                    algorithm = "arpack")
  eigen <- eigen$vectors[,seq_len(n.dim)]
  rownames(eigen) <- barcodes
  colnames(eigen) <- paste0("Trex_", seq_len(ncol(eigen)))
  return(eigen)
}


#Returns appropriate model for autoencoder
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

#Selects columns to normalize input data based on the inputs to the model
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


#Add the eigenvectors to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
adding.DR <- function(sc, reduction, reduction.name) {
  #clusters <- NULL
  #if(length(reduction) > 1) {
  #  clusters <- reduction[[2]]
  #  reduction <- reduction[[1]]
 # }
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
  #No built in clustering for now
  #if(!is.null(clusters)) {
   # clusters <- as.data.frame(clusters)
  #  colnames(clusters) <- paste0(reduction.name, ".cluster")
   # rownames(clusters) <- rownames(grabMeta(sc))
   # sc <- add.meta.data(sc, clusters, colnames(clusters))
  #}
  return(sc)
  
}

AF.col <- c(2,3,4,5,6)
KF.col <- c(7,8,9,10,11,12,13,14,15,16)

#Generates the 30 vector based on autoencoder model 
#First normalizes the value by the min and max of the autoencoder training data
auto.embedder <- function(array.reshape, aa.model, local.max, local.min) {
  for(i in seq_len(ncol(array.reshape))) {
    (array.reshape[,i] - local.min[i])/(local.max[i] - local.min[i])
  }
  array.reshape[is.na(array.reshape)] <- 0
  score <- stats::predict(aa.model, array.reshape)
  return(score)
}

#Code from https://stackoverflow.com/questions/57282842/how-to-efficiently-extract-a-row-or-column-from-a-dist-distance-matrix?rq=1
f <- function (i, j, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
  k <- (2 * n - j) * (j - 1) / 2 + (i - j)
  k[!valid] <- NA_real_
  k
}

#Code from https://stackoverflow.com/questions/57282842/how-to-efficiently-extract-a-row-or-column-from-a-dist-distance-matrix?rq=1
SliceExtract_dist <- function (dist_obj, k) {
  if (length(k) > 1) stop("The function is not 'vectorized'!")
  n <- attr(dist_obj, "Size")
  if (k < 1 || k > n) stop("k out of bound!")
  ##
  i <- 1:(k - 1)
  j <- rep.int(k, k - 1)
  v1 <- dist_obj[f(j, i, dist_obj)]
  ## 
  i <- (k + 1):n
  j <- rep.int(k, n - k)
  v2 <- dist_obj[f(i, j, dist_obj)]
  ## 
  c(v1, 0, v2)
}

neighbor.manager <- function(.row, metric, .length, .j, .nearest.method, .near.neighbor, .threshold, .clone.proportion, .TR) {
  if (metric == "distance") {
      for (k in seq_len(length(.row))) {
        suppressWarnings(.row[k] <- 1- (.row[k]/(.length[.j] + .length[k])/2))
      }
  } else if (metric == "aa.property") {
      max <- max(.row, na.rm = TRUE)
      .row <- (max-abs(.row))/max
  }
  if (.nearest.method == "threshold") {
      neighbor <- which(.row >= .threshold)
  } else if (.nearest.method == "nn") {
      neighbor <- order(.row, decreasing = TRUE)[seq_len(.near.neighbor)]
      neigh.check <- which(.row > 0.99)  
      if (length(neigh.check) > .near.neighbor) {
          matches <- sample(neigh.check, round(.near.neighbor*.clone.proportion))
          .row[.row > 0.99] <- 0
          #Generate "neighborhood"
          close.matches <- order(.row, decreasing = TRUE)
          unique.matches <- unique(.TR[close.matches])[seq_len(.near.neighbor*(1-.clone.proportion))]
          close.matches <- lapply(unique.matches, FUN = function(x) {
              x <- which(.TR == x)
              x <- x[sample(length(x), 1)] })
          neighbor <- c(matches, unlist(close.matches))
      }
    }
  return(neighbor)
}



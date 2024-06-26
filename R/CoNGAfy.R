#' Reduce a single-cell object to representative cells
#' 
#' Generate a single-cell object that has a representation of
#' RNA expression by clonotype. This approach was first introduced
#' in \href{https://pubmed.ncbi.nlm.nih.gov/34426704/}{CoNGA} and was 
#' adapted to R. Please read and cite the author's work. 
#' 
#' @examples
#' trex_values <- CoNGAfy(trex_example, 
#'                        method = "dist",
#'                        features = NULL)
#'                         
#' @param sc Single Cell Object in Seurat or Single Cell Experiment format
#' @param method "mean" or "dist" Generate a mean value across features by clonotype or 
#' use the PCA reduction to identify the cell with the minimal euclidean distance from the
#' clonotype group.
#' @param features Selected genes for the reduction DEFAULT: null will return all genes
#' @param assay The name of the assay or assays to return.
#' @param meta.carry Variables to carry over from the meta data of the single-cell object
#' 
#' @export
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom SummarizedExperiment assay assay<- colData<-
#' 
#' @return Single-cell Object with one cell representing a single clone
#' 
#' 
CoNGAfy <- function(sc, 
                    method = "dist", 
                    features = NULL, 
                    assay = "RNA", 
                    meta.carry = c("CTaa", "CTgene")) {
    sc <- filter.object(sc)
    conga <- NULL
    if(method == "mean") {
        for (x in seq_along(assay)) {
            conga[[x]] <- CoNGA.mean(sc, features, assay[x])
            
        }
    } else if(method == "dist") {
        for (x in seq_along(assay)) {
            conga[[x]] <- CoNGA.dist(sc, features, assay[x])
            
        }
        
    }
    names(conga) <- assay
    if (inherits(x=sc, what ="Seurat")) {
        sc.output <- CreateSeuratObject(conga[[1]], assay = names(conga)[1], project = "Trex")
        if(length(conga) > 1) {
            for(y in 2:length(conga)) {
                sc.output[[names(conga)[y]]] <- CreateAssayObject(conga[[y]])
            }
        }
        CTge <- unique(sc[[]][,c(meta.carry)])
    } else if (inherits(x=sc, what ="SingleCellExperiment")) {
        sc.output <- SingleCellExperiment(assay = conga[[1]])
        if(length(conga) > 1) {
            for(y in 2:length(conga)) {
                assay(sc.output, names(conga)[y]) <- conga[[y]]
            }
        }
        sc.output$CTaa <- rownames(sc.output@colData)
        CTge <- data.frame(unique(sc@colData[,c(meta.carry)]))
    }
    CTge <- meta.handler(CTge, meta.carry)
    clones <- unique(CTge$CTaa)
    rownames(CTge) <- clones
    colnames(CTge) <- meta.carry
    sc.output <- add.meta.data(sc.output, CTge, colnames(CTge))
    return(sc.output)
}



#For all single clones, will use true RNA scaled values
#For multiplets will use the cell with the minimal distance in PCA
#For doublets, will automatically select the first cell. 
#' @importFrom SummarizedExperiment assay
CoNGA.dist <- function(sc, features, assay) {
    if (inherits(x=sc, what ="Seurat")) {
        if(assay == "RNA") {
            data.use <- sc[["pca"]]@cell.embeddings
        } else if (assay == "ADT") {
            data.use <- sc[["apca"]]@cell.embeddings
        }
    } else if (inherits(x=sc, what ="SingleCellExperiment")){
        data.use <- reducedDim(sc, "PCA")
    }
    meta <- grabMeta(sc)
    data <- as.data.frame(meta[,"CTaa"])
    colnames(data) <- "CTaa"
    rownames(data) <- rownames(meta)
    all.clones <- table(data$CTaa)
    unique.clones <- all.clones[which(all.clones > 1)]
    single.clones <- all.clones[all.clones %!in% unique.clones]
    barcodes <- rownames(data)[which(data$CTaa %in% names(single.clones))]
    for (i in seq_along(unique.clones)) {
        loc <- which(data$CTaa == names(unique.clones)[i])
        dist <- as.matrix(dist(data.use[loc,]))
        cell <- names(which(rowSums(dist) == min(rowSums(dist))))
        if (length(cell) > 1) {
            cell <- cell[1]
        }
        barcodes <- c(barcodes, cell)
    }
    if (inherits(x=sc, what ="Seurat")) {
        assay.use <- sc[[assay]]@data
    } else if (inherits(x=sc, what ="SingleCellExperiment")){
        assay.use <- assay(sc)
    }
    features.to.avg <- features %||% rownames(x = assay.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = assay.use))
    data.return <- assay.use[rownames(assay.use) %in% features.assay, colnames(assay.use) %in% barcodes]
    colnames(data.return) <- data$CTaa[match(barcodes, rownames(data))]
    return(data.return)
}
# Adapted from the AverageExpression() function in Seurat
#' @importFrom rlang %||%
#' @importFrom Matrix sparse.model.matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom stats as.formula
CoNGA.mean <- function(sc, features, assay) {
    
    if (inherits(x=sc, what ="Seurat")) {
        data.use <- sc[[assay]]@data
        data.use <- expm1(x = data.use)
    } else if (inherits(x=sc, what ="SingleCellExperiment")){
        data.use <- assay(sc, name = assay)
    }
    features.to.avg <- features %||% rownames(x = data.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
    meta <- grabMeta(sc)
    data <- as.data.frame(meta[,"CTaa"])
    colnames(data) <- "CTaa"
    rownames(data) <- rownames(meta)
    
    #Modified to drop any cells without CTaa
    data.use <- data.use[,!is.na(data[,"CTaa"])]
    data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = FALSE]
    
    for (i in seq_len(ncol(x = data))) {
        data[, i] <- as.factor(x = data[, i])
    }
    num.levels <- sapply(X = seq_len(ncol(x = data)), FUN = function(i) { 
        length(x = levels(x = data[, i]))
    })
    category.matrix <- sparse.model.matrix(object = as.formula(
        object = paste0(
            '~0+',
            paste0(
                "data[,",
                1:length(x = "CTaa"),
                "]",
                collapse = ":"
            )
        )
    ))
    colsums <-Matrix::colSums(x = category.matrix)
    category.matrix <- category.matrix[, colsums > 0]
    colsums <- colsums[colsums > 0]
    
    category.matrix <- sweep(
        x = category.matrix,
        MARGIN = 2,
        STATS = colsums,
        FUN = "/")
    colnames(x = category.matrix) <- sapply(
        X = colnames(x = category.matrix),
        FUN = function(name) {
            name <- gsub(pattern = "data\\[, [1-9]*\\]", replacement = "", x = name)
            return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))), collapse = "_"))
        })
    data.return <- data.use %*% category.matrix
    return(data.return)
}

#' @importFrom stringr str_sort
meta.handler <- function(meta, meta.carry) {
  unique.clones <- unique(meta[,"CTaa"])
  duplicated.clones <- na.omit(unique(meta[,"CTaa"][which(duplicated(meta[,"CTaa"]))]))
  new.meta <- NULL
  for (i in seq_along(duplicated.clones)) {
    meta.tmp <- meta[meta[,"CTaa"] == duplicated.clones[i], ]
    concat.strings <- lapply(meta.carry, function(x) {
            paste0(str_sort(unique(na.omit(meta.tmp[,x]))), collapse = ";")
    })
    new.meta <- rbind(new.meta, unlist(concat.strings))
  }
  old.meta <- meta[meta[,"CTaa"] %!in% duplicated.clones, meta.carry]
  if (length(new.meta) != 0) {
    colnames(new.meta) <- meta.carry
    total.meta <- rbind(old.meta, new.meta)
  } else {
    total.meta <- old.meta
  }
  return(total.meta)
}

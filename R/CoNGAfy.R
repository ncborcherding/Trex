#' Reduce a single-cell object to representative cells
#' 
#' Generate a single-cell object that has a representation of
#' gene expression by clonotype. This approach was first introduced
#' in \href{https://pubmed.ncbi.nlm.nih.gov/34426704/}{CoNGA} and was 
#' adapted to R. Please read and cite the authors work. 
#' 
#' @examples
#' trex_values <- CoNGAfy(trex_example, 
#'                         method = "mean",
#'                         features = NULL)
#'                         
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param method "mean" or "dist" Generate a mean value across features by clonotype or 
#' use the PCA rduction to identify the cell with the minimal euclidena distance from the
#' clonotype group.
#' @param features Selected genes for the reduction DEFAULT: null will return all genes
#' 
#' @export
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @return TSingle-cell Object with 1 cell representing 1 clone
#' 
#' 
CoNGAfy <- function(sc, 
                    method = "mean", 
                    features = NULL) {
    if(method == "mean") {
        conga <- CoNGA.mean(sc, features)
    } else if(methods == "dist") {
        conga <- CoNGA.dist(sc, features)
    }
    if (inherits(x=sc, what ="Seurat")) {
        sc.output <- CreateSeuratObject(conga, project = "Trex")
        CTge <- unique(sc[[]][,c("CTaa", "CTgene")])
        
    } else if (inherits(x=sc, what ="SingleCellExperiment")) {
        sc.output <- SingleCellExperiment(conga)
        sc.output$CTaa <- rownames(sc.output@colData)
        CTge <- data.frame(unique(sc@colData[,c("CTaa", "CTgene")]))
        
    }
    CTge <- CTge[!duplicated(CTge$CTaa),]
    clones <- unique(CTge$CTaa)
    rownames(CTge) <- clones
    colnames(CTge) <- c("CTaa", "CTgene")
    sc.output <- add.meta.data(sc.output, CTge, "CTgene")
    return(sc.output)
}

#For all single clones, will use true RNA scaled values
#For multiplets will use the cell with the minimal distance in PCA
#For doublets, will automatically select the first cell. 
CoNGA.dist <- function(sc, features) {
    if (inherits(x=sc, what ="Seurat")) {
        data.use <- sc[["pca"]]@cell.embeddings
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
        RNA.use <- sc[["RNA"]]@data
    } else if (inherits(x=sc, what ="SingleCellExperiment")){
        RNA.use <- assay(sc)
    }
    features.to.avg <- features %||% rownames(x = RNA.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = RNA.use))
    data.return <- RNA.use[rownames(RNA.use) %in% features.assay, colnames(RNA.use) %in% barcodes]
    colnames(data.return) <- data$CTaa[match(barcodes, rownames(data))]
    return(data.return)
}

#' @importFrom rlang %*% %||%
#' @importFrom Matrix sparse.model.matrix
#' Adapted from the AverageExperssion() function in Seurat
CoNGA.mean <- function(sc, features) {
    if (inherits(x=sc, what ="Seurat")) {
        data.use <- sc[["RNA"]]@data
    } else if (inherits(x=sc, what ="SingleCellExperiment")){
        data.use <- assay(sc)
    }
    
    features.to.avg <- features %||% rownames(x = data.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
    meta <- grabMeta(sc)
    data <- as.data.frame(meta[,"CTaa"])
    colnames(data) <- "CTaa"
    rownames(data) <- rownames(meta)
    data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
    for (i in 1:ncol(x = data)) {
        data[, i] <- as.factor(x = data[, i])
    }
    num.levels <- sapply(X = 1:ncol(x = data), FUN = function(i) { 
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
    colsums <- colSums(x = category.matrix)
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
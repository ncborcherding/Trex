#' Main Trex interface
#' 
#' Use this to run the Trex algorithm to return latent vectors in 
#' the form of a matrix or if you prefer, a maTrex.
#' 
#' @examples
#' trex_values <- maTrex(trex_example, 
#'                         chains = "both", 
#'                         edit.method = "lv",
#'                         AA.properties = "AF",
#'                         AA.method = "auto",
#'                         n.trim = 0,
#'                         c.trim = 0,
#'                         nearest.method = "nn",
#'                         threshold = 0.85,
#'                         near.neighbor = 40,
#'                         add.INKT = FALSE,
#'                         add.MAIT = FALSE, 
#'                         n.dim = 40,
#'                         species = "human")
#'                         
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{stringdist} or NULL
#' to not include edit distance layer in multiplex network
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or NULL to not 
#' include edit distance layer in multiplex network
#' @param AA.method method for deriving pairwise distances from amino acid properties, 
#' either "mean" or the average amino acid characteristics for given cdr3 sequence or "auto"
#' for using an autoencoder model. The autoencoder selection  will not allow for trimming of 
#' the cdr3 sequences.
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' @param nearest.method method for defining neighbors, either "threshold" for normalized
#' distances above a certain value or "nn" for nearest neighbor by normalized distances. 
#' @param threshold If nearest.method = "threshold" -the value for distance measures 
#' above which simplify for adjacency matrix
#' @param near.neighbor If nearest.method = "nn" - the number of nearest neighbors to use 
#' in generating an adjacency matrix. If the number of copies of a clone is greater than the
#' near.neighbor value used, the near.neighbor will randomly sample unique clones that are nearest
#' to the clone being compared.
#' @param clone.proportion The proportion of nearest neighbors to return from the a single clone
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param n.dim The number of Trex dimensions to return, similar to PCA dimensions
#' @param species Indicate "human" or "mouse" for gene-based metrics
#' @param seed seed for the random number generator
#' 
#' @export
#' @importFrom SeuratObject CreateDimReducObject
#' 
#' @return Trex eigenvectors calculated from multiplex network
maTrex <- function(sc, 
                    chains = "both", 
                    edit.method = "lv",
                    AA.properties = "AF",
                    AA.method = "auto",
                    n.trim = 0,
                    c.trim = 0,
                    nearest.method = "nn",
                    threshold = 0.85,
                    near.neighbor = 40,
                    clone.proportion = 0.5,
                    add.INKT = TRUE,
                    add.MAIT = TRUE, 
                    n.dim = 40,
                    species = "human", 
                    seed = 42) {
    set.seed(seed)
    TCR <- getTCR(sc, chains)
    print("Calculating the Edit Distance for CDR3 AA sequence...")
    if (!is.null(edit.method)) {
        network <- distanceMatrix(TCR, edit.method, nearest.method, near.neighbor, threshold, c.trim, n.trim, clone.proportion, return.dims = FALSE)
    } else {
        network <- NULL
    }
    
    if ((AA.properties %in% c("AF", "KF", "both", "all"))[1]) {
        print("Calculating the Amino Acid Properties...")
        AA.knn <- aaProperty(TCR, c.trim, n.trim, nearest.method, near.neighbor, threshold, AA.method, AA.properties, clone.proportion, return.dims = FALSE)
        network <- c(network, AA.knn)
    }
    
    if (add.INKT) {
        print("Calculating the INKT gene usage...")
        tmpscore <- scoreINKT(TCR, species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            network <- c(network, tmp.knn)
            names(network)[length(network)] <- "INKT"
        }
    }
    if (add.MAIT) {
        print("Calculating the MAIT gene usage...")
        tmpscore <- scoreMAIT(TCR, species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            network <- c(network, tmp.knn)
            names(network)[length(network)] <- "MAIT"
        }
    }
    print("Calculating Latent Vectors from multiplex network...")
    barcodes <- rownames(grabMeta(sc))
    reduction <- multiplex.network(network, n.dim, barcodes)
    return(reduction)
}

#' Trex single cell calculation
#'
#'Run Trex algorithm with Seurat or SingleCellExperiment pipelines
#'
#' @examples
#' trex_example <- runTrex(trex_example, 
#'                         AA.properties = "AF", 
#'                         nearest.method = "nn",
#'                         near.neighbor = 40,
#'                         threshold = 0.85,
#'                         reduction.name = "Trex.AF")
#'                         
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{stringdist}
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' "AF" = Atchley factors, "KF" = Kidera factors, or "both"
#' @param AA.method method for deriving pairwise distances from amino acid properties, 
#' either "mean" or the average amino acid characteristics for given cdr3 sequence or "auto"
#' for using an autoencoder model. The autoencoder selection  will not allow for trimming of 
#' the cdr3 sequences.
#' @param reduction.name Keyword to save Trex reduction. Useful if you want
#' to try Trex with multiple parameters 
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' @param nearest.method method for defining neighbors, either "threshold" for normalized
#' distances above a certain value or "nn" for nearest neighbor by normalized distances. 
#' @param threshold If nearest.method = "threshold" -the value for distance measures 
#' above which simplify for adjacency matrix
#' @param near.neighbor If nearest.method = "nn" - the number of nearest neighbors to use 
#' in generating an adjacency matrix. If the number of copies of a clone is greater than the
#' near.neighbor value used, the near.neighbor will randomly sample unique clones that are nearest
#' to the clone being compared.
#' @param clone.proportion The proportion of nearest neighbors to return from the a single clone
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param n.dim The number of Trex dimensions to return, similar to PCA dimensions
#' @param species Indicate "human" or "mouse" for gene-based metrics
#' @param seed seed for the random number generator
#' @export
#' @return Seurat or SingleCellExperiment object with Trex dimensions placed 
#' into the dimensional reduction slot. 
#' 
runTrex <- function(sc, 
                    chains = "both", 
                    edit.method = "lv",
                    AA.properties = "AF",
                    AA.method = "auto",
                    n.trim = 0,
                    c.trim = 0,
                    nearest.method = "nn",
                    threshold = 0.85,
                    near.neighbor = 40,
                    clone.proportion = 0.5,
                    add.INKT = TRUE,
                    add.MAIT = TRUE, 
                    n.dim = 40,
                    species = "human",
                    reduction.name = "Trex",
                    seed = 42) {

    cells.chains <- rownames(sc[[]][!is.na(sc[["cloneType"]]),])
    sc <- subset(sc, cells = cells.chains)
    reduction <- maTrex(sc,
                        chains, 
                        edit.method,
                        AA.properties,
                        AA.method,
                        n.trim, 
                        c.trim,
                        nearest.method,
                        threshold,
                        near.neighbor,
                        clone.proportion,
                        add.INKT,
                        add.MAIT,
                        n.dim,
                        species, 
                        seed)
    TCR <- getTCR(sc, chains)
    if (add.INKT) {
        tmpscore <- scoreINKT(TCR, species)
        sc <- add.meta.data(sc, tmpscore, "IKNT.score")
    }
    if (add.MAIT) {
        tmpscore <- scoreMAIT(TCR, species)
        sc <- add.meta.data(sc, tmpscore, "MAIT.score")
    }
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}

#' Remove TCR genes from variable gene results
#'
#'Most single-cell workflows use highly-expressed and highly-variable
#'genes for the initial calculation of PCA and subsequent dimensional
#'reduction. This function will remove the TCR genes from the variable
#'features in the Seurat object or from a vector generated by
#'the Bioconductor scran workflow. 
#'
#' @examples
#' x <- trex_example
#' x <- quietTCRgenes(x)
#' 
#' @param sc Single Cell Object in Seurat format or vector of variable genes to use in reduction
#' @export
#' @return Seurat object or vector list with TCR genes removed.
quietTCRgenes <- function(sc) {
    unwanted_genes <- "TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*"
    if (inherits(x=sc, what ="Seurat")) {
        unwanted_genes <- grep(pattern = unwanted_genes, x = sc[["RNA"]]@var.features, value = TRUE)
        sc[["RNA"]]@var.features <- sc[["RNA"]]@var.features[sc[["RNA"]]@var.features %!in% unwanted_genes]
    } else {
        #Bioconductor scran pipelines uses vector of variable genes for DR
        unwanted_genes <- grep(pattern = unwanted_genes, x = sc, value = TRUE)
        sc <- sc[sc %!in% unwanted_genes]
    }
    return(sc)
}

#' Cluster clones using the Trex dimensional reductions
#' 
#' Use this to return clusters for clonotypes based on 
#' the \link[bluster]{bluster} clustering paramaters.
#' 
#' @examples
#' \dontrun{
#' sc <- clonalCommunity(sc, 
#'                       reduction.name = NULL, 
#'                       cluster.parameter = NNGraphParam())
#' }
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param reduction.name Name of the dimensional reduction output from runTrex()
#' @param cluster.parameter The community detection algorithm in \link[bluster]{bluster}
#' @param ... For the generic, further arguments to pass to specific methods.
#' @importFrom bluster clusterRows NNGraphParam HclustParam KmeansParam KNNGraphParam PamParam SNNGraphParam SomParam
#' @export
#' @return Single-Cell Object with trex.clusters in the meta.data
clonalCommunity <- function(sc, 
                            reduction.name = NULL, 
                            cluster.parameter=NNGraphParam(k=30, ...), 
                            ...) {
    if (inherits(x=sc, what ="Seurat")) { 
        dim.red <- sc[[reduction.name]] 
        dim.red <- dim.red@cell.embeddings
    } else {
        dim.red <- reducedDim(sc, reduction.name)
    }
    clusters <- clusterRows(dim.red, BLUSPARAM=cluster.parameter)
    clus.df <- data.frame("trex.clusters" = paste0("trex.", clusters))
    rownames(clus.df) <- rownames(dim.red)
    sc <- add.meta.data(sc, clus.df, colnames(clus.df))
}

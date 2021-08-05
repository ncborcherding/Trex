#' Main Trex interface
#' 
#' Use this to run the Trex algorithm to return latent vectors

#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{sstringdist} or NULL
#' to not include edit distance layer in multiplex network
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or NULL to n
#' ot include edit distance layer in multiplex network
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
#' near.neighbor value used, the near.neighbor will randomly return neighbors/clones. 
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param n.dim The number of Trex dimensions to return, similar to PCA dimensions
#' @param species Indicate "human" or "mouse" for gene-based metrics
#' 
#' @importFrom SeuratObject CreateDimReducObject
#' 
#' @return Trex eigen vectors caculated from multiplex network
maTrex <- function(sc, 
                    chains = "both", 
                    edit.method = "lv",
                    AA.properties = c("AF", "KF", "both"),
                    AA.method = "auto",
                    n.trim = 0,
                    c.trim = 0,
                    nearest.method = "threshold",
                    threshold = 0.85,
                    near.neighbor = NULL,
                    add.INKT = TRUE,
                    add.MAIT = TRUE, 
                    n.dim = 30,
                    species = "human") {
    TCR <- getTCR(sc, chains)
    print("Calculating the Edit Distance for CDR3 AA sequence...")
    if (!is.null(edit.emthod)) {
        network <- distanceMatrix(TCR, edit.method, nearest.method, threshold, near.neighbor, c.trim, n.trim)
    } else {
        network <- NULL
    }
    
    if (unique(c("AF", "KF", "both", "all") %in% AA.properties)[1]) {
        print("Calculating the Amino Acid Properties...")
        AA.knn <- aaProperty(TCR, c.trim, n.trim, nearest.method, threshold, near.neighbor, AA.method, AA.properties)
        network <- add.to.network(network, AA.knn, paste0(names(TCR), ".AA")) 
    }
    
    multi.network <- list()
    for (i in seq_along(network)) {
        multi.network[[i]] <- get.knn(TCR[[1]]$barcode, out_matrix = network[[i]], 
                                       nearest.method, near.neighbor, threshold)
    }
    names(multi.network) <- names(network)
    
    if (add.INKT) {
        print("Calculating the INKT gene usage...")
        tmpscore <- scoreINKT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
    }
    if (add.MAIT) {
        print("Calculating the MAIT gene usage...")
        tmpscore <- scoreMAIT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
    }
    print("Calculating Latent Vectors from multiplex network...")
    barcodes <- rownames(grabMeta(sc))
    reduction <- multiplex.network(multi.network, n.dim, barcodes)
    return(reduction)
}

#' Trex single cell calculation
#'
#'Run Trex algorithm with Seurat or SingleCellExperiment pipelines


#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{sstringdist}
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' AF = Atchley factors, KF = Kidera factors, or other (polarity, molecular weight, etc), 
#' or none
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
#' near.neighbor value used, the near.neighbor will randomly return neighbors/clones. 
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param n.dim The number of Trex dimensions to return, similar to PCA dimensions
#' @param species Indicate "human" or "mouse" for gene-based metrics

#' @importFrom SeuratObject subset
#' 
#' #' @return Seurat or SingleCellExperiment object with Trex dimensions placed 
#' into the dimensional reduction slot. 
#' 
runTrex <- function(sc, 
                   chains = "both", 
                   edit.method = "lv",
                   AA.properties = c("AF", "KF", "other"),
                   AA.method = "auto",
                   reduction.name = "Trex",
                   n.trim = 0,
                   c.trim = 0,
                   nearest.method = "threshold",
                   threshold = 0.85,
                   near.neighbor = NULL,
                   add.INKT = TRUE,
                   add.MAIT = TRUE, 
                   n.dim = 30,
                   species = "human") {
        
    cells.chains <- rownames(sc[[]][!is.na(sc[["cloneType"]]),])
    sc <- subset(sc, cells = cells.chains)
    reduction <- maTrex(sc,
                        chains, 
                        edit.method,
                        AA.properties,
                        AA.method,
                        c.trim,
                        n.trim,
                        nearest.method,
                        threshold,
                        near.neighbor,
                        add.INKT,
                        add.MAIT,
                        n.dim,
                        species)
    TCR <- getTCR(sc, chains)
    if (add.INKT) {
        tmpscore <- scoreINKT(TCR, species = species)
        sc <- add.meta.data(sc, tmpscore, "IKNT.score")
    }
    if (add.MAIT) {
        tmpscore <- scoreMAIT(TCR, species = species)
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

#' 
#' @param sc Single Cell Object in Seurat format or vector of variable genes to use in reduction
#' @export
#' @return Seurat object or vector list with TCR genes removed.
quietTCRgenes <- function(sc) {
    unwanted_genes <- "TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*"
    if (inherits(x=sc, what ="Seurat")) {
        unwanted_genes <- grep(pattern = unwanted_genes, x = sc[["RNA"]]@var.features, value = T)
        sc[["RNA"]]@var.features <- sc[["RNA"]]@var.features[sc[["RNA"]]@var.features %!in% unwanted_genes]
    } else {
        #Biocondutor scran pipelines uses vector of variable genes for DR
        unwanted_genes <- grep(pattern = unwanted_genes, x = sc, value = T)
        sc <- sc[sc %!in% unwanted_genes]
    }
    return(sc)
}


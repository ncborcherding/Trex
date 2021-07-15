#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{sstringdist}
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' AF = Atchley factors, KF = Kidera factors, or other (polarity, molecular weight, etc), 
#' or none
#' @param AA.method method for deriving pairwise distances from amino acid properties, 
#' either "mean" or the average amino acid characteristics for given cdr3 sequence or "auto"
#' for using an auto-encoder model.
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param nearest.method method for defining neighbors, either "threshold" for normalized
#' distances above a certain value or "nn" for nearest neighbor by normalized distances. 
#' @param threshold If nearest.method = "threshold" -the value for distance measures 
#' above which simplify for adjacency matrix
#' @param near.neighbor If nearest.method = "nn" - the number of nearest neighbors to use 
#' in generating an adjacency matrix. If the number of copies of a clone is greater than the
#' near.neighbor value used, the near.neighbor will randomly return neighbors/clones. 
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param species Indicate "human" or "mouse" for gene-based metrics
#' 
#How to handle non-T-cell
#CD4 to CD8 model
#AA - multiplex factors  vs total
#   - does the threshold make sense? 
# KNN vs threshold neighbor
#' @importFrom SeuratObject CreateDimReducObject

runTrex <- function(sc, 
                   chains = "both", 
                   edit.method = "lv",
                   AA.properties = c("AF", "KF", "other"),
                   AA.method = "auto",
                   c.trim = 0,
                   n.trim = 0,
                   nearest.method = "threshold",
                   threshold = 0.85,
                   near.neighbor = NULL,
                   add.INKT = TRUE,
                   add.MAIT = TRUE, 
                   species = "human") {
    #Add Filtering step
    TCR <- getTCR(sc, chains)
    print("Calculating the Edit Distance for CDR3 AA sequence...")
    multi.network <- distanceMatrix(TCR, edit.method, c.trim, n.trim, nearest.method, threshold, near.neighbor)
    
    if (unique(c("AF", "KF", "other") %in% AA.properties)[1]) {
        print("Calculating the Amino Acid Properties...")
        if(AA.method == "mean") {
            AA.knn <- aaProperty(TCR, c.trim, n.trim, AA.properties)
        } else if (AA.method == "auto") {
            AA.knn <- aaAutoEncoder(TCR, c.trim, n.trim, AA.properties)
        }
        multi.network <- add.to.network(multi.network, AA.knn, paste0(names(TCR), ".AA")) 
    }
    if (add.INKT) {
        print("Calculating the INKT gene usage...")
        tmpscore <- scoreINKT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
            barcodes <- tmpscore$barcode
            tmpscore <- as.data.frame(tmpscore[,2])
            colnames(tmpscore) <- "IKNT.score"
            rownames(tmpscore) <- barcodes
            sc <- add.meta.data(sc, tmpscore)
        }
    }
    if (add.MAIT) {
        print("Calculating the MAIT gene usage...")
        tmpscore <- scoreMAIT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
            barcodes <- tmpscore$barcode
            tmpscore <- as.data.frame(tmpscore[,2])
            colnames(tmpscore) <- "MAIT.score"
            rownames(tmpscore) <- barcodes
            sc <- add.meta.data(sc, tmpscore)
        }
    }
    print("Multiplexing Nodes into single graph...")
    graph <- multiplex.network(multi.network)
    return(graph)
}

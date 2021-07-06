#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA, TRB, both
#' @param edit.method distance measures inherited from \link[stringdist]{sstringdist}
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' AF = Atchley factors, KF = Kidera factors, or other (polarity, molecular weight, etc), 
#' or none
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param threshold The value for distance measures above which simplify for adjacency matrix
#' @param add.INKT Add a additional layer for invariant natural killer T cells based on genes
#' @param add.MAIT Add a additional layer for Mucosal-associated invariant T cells based on genes
#' @param species Indicate "human" or "mouse" for gene-based metrics
#' 
#How to handle non-T-cell
#CD4 to CD8 model
#AA - multiplex factors  vs total
#   - pearson correlation used for distances
#   - does the threshold make sense? 
# KNN vs threshold neighbor
#' @importFrom SeuratObject CreateDimReducObject

runTrex <- function(sc, 
                   chains = "both", 
                   edit.method = "lv",
                   AA.properties = c("AF", "KF", "other"),
                   c.trim = 0,
                   n.trim = 0,
                   threshold = 0.85,
                   add.INKT = TRUE,
                   add.MAIT = TRUE, 
                   species = "human") {
    TCR <- getTCR(sc, chains)
    print("Calculating the Edit Distance for CDR3 AA sequence...")
    multi.network <- distanceMatrix(TCR, edit.method, c.trim, n.trim, threshold)
    
    if (AA.properties %in% c("AF", "KF", "other")) {
        print("Calculating the Amino Acid Properties...")
        AA.knn <- aaProperty(TCR, c.trim, n.trim, threshold, AA.propertiess)
        multi.network <- add.to.network(multi.network, AA.knn, names(TCR)) 
    }
    if (add.INKT) {
        print("Calculating the INKT gene usage...")
        tmpscore <- scoreINKT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
        barcodes <- tmpscore$barcode
        tmpscore <- as.data.frame(tmpscore[,2])
        colnames(tmpscore) <- "IKNT.score"
        rownames(tmpscore) <- barcodes
        sc <- add.meta.data(sc, tmpscore)
    }
    if (add.MAIT) {
        print("Calculating the MAIT gene usage...")
        tmpscore <- scoreMAIT(TCR, species = species)
        if (length(which(tmpscore$score > 0)) != 0) {
            tmp.knn <- gene.to.knn(tmpscore)
            multi.network <- add.to.network(multi.network, tmp.knn, "INKT") 
        }
        barcodes <- tmpscore$barcode
        tmpscore <- as.data.frame(tmpscore[,2])
        colnames(tmpscore) <- "MAIT.score"
        rownames(tmpscore) <- barcodes
        sc <- add.meta.data(sc, tmpscore)
    }
    print("Multiplexing Nodes into single graph...")
    adj.matrix <- multiplex.network(multi.network)
}
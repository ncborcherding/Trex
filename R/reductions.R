#' Reduce TCR data based on amino acid characteristics
#' 
#' Use this to get the vector returned from the TCR autoencoder 
#' or mean values based on indicated amino acid property before
#' the network generation.
#' 
#' @examples
#' AA.dim  <- aaReduction(trex_example, 
#'                         chains = "TRB", 
#'                         AA.properties = "AF",
#'                         AA.method = "auto",
#'                         n.trim = 0,
#'                         c.trim = 0)
#' 
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA or TRB (no both option)
#' @param AA.properties Amino acid properties to use for distance calculation: 
#' "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or NULL to not 
#' include edit distance layer in multiplex network
#' @param AA.method method for deriving pairwise distances from amino acid properties, 
#' either "mean" or the average amino acid characteristics for given cdr3 sequence or "auto"
#' for using an autoencoder model. The autoencoder selection  will not allow for trimming of 
#' the cdr3 sequences.
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' 
#' @export
#' 
#' @return 30-dimensional autenocoder bottleneck layer or mean values for each cell
aaReduction <- function(sc, 
                        chains = "TRB", 
                        AA.properties = "AF",
                        AA.method = "auto",
                        n.trim = 0,
                        c.trim = 0) {
    if (chains %!in% c("TRB", "TRA")) {
        stop("Please Select either TRA or TRB as the chain to return")
    }
    TCR <- getTCR(sc, chains)
    AA <- aaProperty(TCR, 
                         n.trim, 
                         c.trim, 
                         nearest.method = "nn",
                         threshold = 0.85, 
                         near.neighbor = 40, 
                         AA.method, 
                         AA.properties, 
                         return.dims = TRUE)
    return(AA)
}

#' Reduce TCR data based on pairwise edit distance
#' 
#' Use this to get the pairwise edit distance of the TCRs. Distance methods 
#' available are from \link[stringdist]{stringdist} package.s
#' 
#' @examples
#' dist.dim  <- distReduction(trex_example, 
#'                         chains = "TRB", 
#'                         edit.method = "lv",
#'                         n.trim = 0,
#'                         c.trim = 0)
#' 
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA or TRB (no both option)
#' @param edit.method distance measures inherited from \link[stringdist]{stringdist} or NULL
#' to not include edit distance layer in multiplex network
#' @param n.trim length from N-terminus to trim on the CDR3 aa sequence
#' @param c.trim length from C-terminus to trim on the CDR3 aa sequence
#' 
#' @export
#' 
#' @return distance object based on the indicated index measure


distReduction <- function(sc, 
                        chains = "TRB", 
                        edit.method = "lv",
                        n.trim = 0,
                        c.trim = 0) {
    if (chains %!in% c("TRB", "TRA")) {
        stop("Please Select either TRA or TRB as the chain to return")
    }
    TCR <- getTCR(sc, chains)
    DIST<- distanceMatrix(TCR, 
                          edit.method = "lv",
                          nearest.method = "nn",
                          near.neighbor = 40, 
                          threshold = 0.85, 
                          n.trim, 
                          c.trim, 
                          return.dims = TRUE)
        
        
    return(DIST)
}
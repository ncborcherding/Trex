#' Main Trex interface
#' 
#' Use this to run the Trex algorithm to return latent vectors in 
#' the form of a matrix or if you prefer.
#' 
#' @examples
#'Trex_values <- maTrex(trex_example, 
#'                      chains = "TRA",
#'                      method = "encoder",
#'                      encoder.model = "VAE",
#'                      encoder.input = "AF")
#'                           
#'Trex_values <- maTrex(trex_example, 
#'                      chains = "TRA",
#'                      method = "geometric",
#'                      theta = pi)
#'                         
#' @param sc Single-cell object or
#' the output of combineTCR() in scRepertoire
#' @param chains TCR chain to use - \strong{TRA} or \strong{TRB}.
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on BLOSUM62 matrix
#' @param encoder.model "AE" = dense autoencoder or "VAE" = variation autoencoder
#' @param encoder.input  Atchley factors (\strong{AF}), Kidera factors ((\strong{AF})), 
#' Atchley and Kidera factors ((\strong{both})), or One-Hot Encoder (\strong{OHE}).
#' @param theta angle to use for geometric transformation
#' 
#' @export
#' @importFrom SeuratObject CreateDimReducObject
#' 
#' @return Trex encoded values from the autoencoder
maTrex <- function(sc, 
                   chains = "TRA", 
                   method = "encoder",
                   encoder.model = "VAE", 
                   encoder.input = "AF",
                   theta = pi) {
    TCR <- getTCR(sc, chains)
    checkLength(TCR[[1]])
    if (method == "encoder" && encoder.input %in% c("AF", "KF", "both", "all", "OHE")) {
        print("Calculating the encoding values...")
        reduction <- .encoder(TCR, encoder.input, encoder.model)
    } else if (method == "geometric") {
        print("Performing geometric transformation...")
        reduction <- .geometric.encoding(TCR, theta)
    }
    return(reduction)
}

#' Trex single cell calculation
#'
#'Run Trex algorithm with Seurat or SingleCellExperiment pipelines
#'
#' @examples
#' trex_example <- runTrex(trex_example, 
#'                         chains = "TRA",
#'                         method = "encoder",
#'                         encoder.model = "VAE",
#'                         encoder.input = "AF")
#'                         
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA or TRB
#' @param method "encoder" = using deep learning autoencoders or 
#' "geometric" = geomteric transformations based on BLOSUM62 matrix
#' @param encoder.model "AE" = dense autoencoder or "VAE" = variation autoencoder
#' @param encoder.input "AF" = Atchley factors, "KF" = Kidera factors, "both" = AF and KF, or "OHE" for
#' One Hot Autoencoder
#' @param theta angle to use for geometric transformation
#' @param reduction.name Keyword to save Trex reduction. Useful if you want
#' to try Trex with multiple parameters 
#' @export
#' @return Seurat or SingleCellExperiment object with Trex dimensions placed 
#' into the dimensional reduction slot. 
#' 
runTrex <- function(sc, 
                    chains = "TRA", 
                    method = "encoder",
                    encoder.model = "VAE", 
                    encoder.input = "AF",
                    theta = pi,
                    reduction.name = "Trex") {
    checkSingleObject(sc)
    cells.chains <- rownames(sc[[]][!is.na(sc[["CTaa"]]),])
    sc <- subset(sc, cells = cells.chains)
    reduction <- maTrex(sc = sc,
                        chains = chains, 
                        method = method,
                        encoder.model = encoder.model, 
                        encoder.input = encoder.input,
                        theta = theta)
    TCR <- getTCR(sc, chains)
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}

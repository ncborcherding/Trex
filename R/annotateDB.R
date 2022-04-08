#' Use reference database to annotate cdr3 sequences
#'
#' Use meta database of cdr3 sequences as a reference to 
#' cdr3 aa sequence of the single cells. 
#'
#' @examples
#' trex_example <- annotateDB(trex_example, 
#'                         chains = "TRB")
#' @param sc Single Cell Object in Seurat or SingleCell Experiment format
#' @param chains TRA or TRB
#' @export
#' @return Seurat or SingleCellExperiment object with epitope information
#' added to the metadata. 

annotateDB <- function(sc,
                       chains = "TRB") {
  TCR <- getTCR(sc, chains)
  relevent.db <- as.data.frame(Trex.database[[chains]])
  TCR <- TCR[[1]]
  cdr3.match <- TCR[TCR$cdr3_aa %in% relevent.db[,1],]
  cdr3.match <- merge(cdr3.match, relevent.db, by.x ="cdr3_aa", by.y ="CDR3", all.x = TRUE)
  barcodes <- cdr3.match$barcode
  cdr3.match <- cdr3.match[,c("Epitope.target", "Epitope.sequence", "Epitope.species", "Tissue", "CellType", "Database")]
  colnames(cdr3.match) <- paste0(chains, "_", colnames(cdr3.match))
  rownames(cdr3.match) <- barcodes
  sc <- add.meta.data(sc, cdr3.match)
}

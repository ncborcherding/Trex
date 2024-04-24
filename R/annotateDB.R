#' Use reference database to annotate cdr3 sequences
#'
#' Use meta database of cdr3 sequences as a reference to 
#' cdr3 aa sequence of the single cells. 
#'
#' @examples
#' trex_example <- annotateDB(trex_example, 
#'                         chains = "TRB")
#' @param sc Single Cell Object in Seurat or Single-Cell Experiment format
#' @param chains TCR chain to use - \strong{TRA} or \strong{TRB}
#' @param edit.distance Number of amino acid residues different than reference
#' using Levenshtein distance. Distances of <= 2 are optimal.
#' @importFrom stringdist stringdist
#' @export
#' @return Seurat or SingleCellExperiment object with epitope information
#' added to the metadata. 

annotateDB <- function(sc,
                       chains = "TRB",
                       edit.distance = 0) {
  TCR <- getTCR(sc, chains)
  relevent.db <- as.data.frame(Trex.database[[chains]])
  TCR <- TCR[[1]]
  cdr3.match <- TCR[TCR$cdr3_aa %in% relevent.db[,1],]
  cdr3.match <- merge(cdr3.match, relevent.db, by.x ="cdr3_aa", by.y ="CDR3", all.x = TRUE)
  barcodes <- cdr3.match$barcode
  cdr3.match <- cdr3.match[,c("Epitope.target", "Epitope.sequence", "Epitope.species", "Tissue", "Cell.type", "Database")]
  colnames(cdr3.match) <- paste0(chains, "_", colnames(cdr3.match))
  rownames(cdr3.match) <- barcodes
  if (edit.distance > 0) {
    annotated.cdr3 <- unique(TCR[TCR$barcode %in% barcodes, 2:5])[,1]
    
    additions <- lapply(annotated.cdr3, function(x) {
      pos <- which(stringdist(x, TCR$cdr3_aa, method = "lv") <= edit.distance & 
                     stringdist(x, TCR$cdr3_aa, method = "lv") > 0)
      pos
    })
    for(i in seq_along(additions)) {
      if(length(additions[[i]]) == 0) {
        next()
      }
      new.annotation <- cdr3.match[rownames(cdr3.match) %in% TCR[TCR$cdr3_aa %in% annotated.cdr3[i],]$barcode[1],]
      tmp <- as.data.frame(t(matrix(nrow = ncol(new.annotation), 
                              ncol = length(TCR[additions[[i]],]$barcode),
                              data = unlist(new.annotation))))
      
      rownames(tmp) <- TCR[additions[[i]],]$barcode
      colnames(tmp) <- colnames(cdr3.match)
      cdr3.match <- rbind(cdr3.match, tmp)
    }
  }
  sc <- add.meta.data(sc, cdr3.match)
}

# Transforming the CDR3 sequence into vectors
#' @importFrom keras load_model_hdf5
#' @importFrom reticulate array_reshape
.encoder <- function(TCR, 
                     encoder.input = encoder.input, 
                     encoder.model = encoder.model) { 
  return <- list() ### Need to add reference data
  reference <- Trex.Data[[1]] #AA properties
  col.ref <- grep(tolower(paste(encoder.input, collapse = "|")), colnames(reference))
  length <- NULL
  if (encoder.input == " both") {
    column.ref <- unique(sort(c(AF.col, KF.col)))
  } else {
    column.ref <- unique(sort(col.ref))
  }
  chain <- names(TCR)
  aa.model <- quiet(aa.model.loader(chain, encoder.input, encoder.model))
  membership <- TCR[[1]]
  names <- membership$barcode
  array.reshape <- NULL
    
  cells <- unique(membership[,"barcode"])
  score <- NULL
  for (n in seq_len(length(cells))) {
      tmp.CDR <- membership[membership$barcode == cells[n],]$cdr3_aa
      refer <- unlist(strsplit(tmp.CDR, ""))
      refer <- c(refer, rep(NA, 60 - length(refer)))
      if(encoder.input == "OHE") {
        int <- one.hot.organizer(refer)
        array.reshape.tmp <- array_reshape(int, 1260)
      }else {
        int <- reference[match(refer, reference$aa),c(1,col.ref)]
        int <- as.matrix(int[,-1])
        array.reshape.tmp <- array_reshape(int, length(col.ref)*60)
      }
      score.tmp <- auto.embedder(array.reshape.tmp, aa.model, encoder.input)
      score <- rbind(score, score.tmp)
    }
    score <- data.frame(unique(membership[,"barcode"]), score)
    barcodes <- score[,1]
    score <- score[,-1]
    rownames(score) <- barcodes
    colnames(score) <- paste0("Trex_", seq_len(ncol(score)))
    return(score)
}




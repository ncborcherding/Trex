# Transforming the CDR3 sequence into vectors
#' @importFrom reticulate array_reshape
.encoder <- function(TCR, 
                     encoder.input = encoder.input, 
                     encoder.model = encoder.model) { 

  #Get model path
  chain <- names(TCR)
  model.path <- aa.model.loader(chain, encoder.input, encoder.model)

  #Get sequences
  sequences <- TCR[[1]]$cdr3_aa
  barcodes <- TCR[[1]]$barcode

  #Encode sequences
  property.set <- switch(encoder.input,
                         "AF" = "atchleyFactors",
                         "KF" = "kideraFactors",
                         "OHE" = "OHE")

  if (property.set == "OHE") {
      encoded.values <- onehotEncoder(sequences, max.length = 60)[[2]]
  } else {
      encoded.values <- propertyEncoder(sequences, max.length = 60, property.set = property.set)[[2]]
  }

  #Get reduction
  reduction <- basiliskRun(
    env = TrexEnv,
    fun = function(mpath, xmat) {
      keras <- reticulate::import("keras", delay_load = FALSE)
      model <- NULL
      pred <- NULL

      tryCatch({
        model <- keras$models$load_model(mpath)
        pred  <- model$predict(xmat)
        as.array(pred) # This will be the return value
      }, finally = {
        rm(pred)
        rm(model)
        # Clear the Keras/TensorFlow backend session
        keras$backend$clear_session()
        gc()
      })
    },
    mpath = model.path,
    xmat  = encoded.values
  )

  reduction <- as.data.frame(reduction)
  rownames(reduction) <- barcodes
  colnames(reduction) <- paste0("Trex_", seq_len(ncol(reduction)))
  return(reduction)
}




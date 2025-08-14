"%!in%" <- Negate("%in%")

aa.eval <- function(x) { x %in% c("AF", "KF", "other")}

# Add to meta data some of the metrics calculated
#' @importFrom rlang %||%
#' @importFrom SingleCellExperiment colData
add.meta.data <- function(sc, meta, header) {
if (inherits(x=sc, what ="Seurat")) { 
  col.name <- names(meta) %||% colnames(meta)
  sc[[col.name]] <- meta
} else {
  rownames <- rownames(colData(sc))
  colData(sc) <- cbind(colData(sc), 
          meta[rownames,])[, union(colnames(colData(sc)),  colnames(meta))]
  rownames(colData(sc)) <- rownames  
}
  return(sc)
}


#This is to grab the metadata from a Seurat or SCE object
#' @importFrom SingleCellExperiment colData 
grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}

#Shhhhhh
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#This is to check the single-cell expression object
checkSingleObject <- function(sc) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
}

#This is to check that all the cdr3 sequences are < 60 residues
checkLength <- function(x) {
  if(any(na.omit(nchar(x[,"cdr3_aa"])) > 60)) {
    stop("Models have been trained on cdr3 sequences 
         less than 60 amino acid residues. Please
         filter the larger sequences before running")
  }
}
#Returns appropriate model for autoencoder
#' @importFrom utils download.file read.csv
#' @importFrom tools R_user_dir
aa.model.loader <- function(chain,
                            encoder.input,
                            encoder.model) {

  ## 1. Expected filename
  file_name <- paste0(chain, "_", encoder.input, "_", encoder.model, ".h5")

  ## 2. Sanity-check against metadata.csv
  meta <- read.csv(
    system.file("extdata", "metadata.csv", package = "Trex"),
    stringsAsFactors = FALSE
  )

  if (!file_name %in% meta[[1]])
    stop("Model '", file_name, "' is not listed in metadata.csv.")

  ## 3. Cache directory
  cache_dir <- tools::R_user_dir("Trex", which = "cache")
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  local_path <- file.path(cache_dir, file_name)

  ## 4. Download if we have never seen this model before
  if (!file.exists(local_path)) {
    message("Downloading model '", file_name, " ...")
    #I am using a placeholder for the Zenodo record ID.
    #The user will need to provide the correct ID.
    base_url <- "https://zenodo.org/record/PLACEHOLDER/files"
    status   <- utils::download.file(
      url      = file.path(base_url, file_name),
      destfile = local_path,
      mode     = "wb",
      quiet    = TRUE
    )
    if (status != 0)
      stop("Download of model '", file_name,
           "' failed (status ", status, ").")
  }

  ## 5. Done return the path for use in basiliskRun()
  normalizePath(local_path, winslash = "/")
}


#Add the eigenvectors to single cell object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
adding.DR <- function(sc, reduction, reduction.name) {
  if (inherits(sc, "Seurat")) {
    DR <- suppressWarnings(CreateDimReducObject(
      embeddings = as.matrix(reduction),
      loadings = as.matrix(reduction),
      projected = as.matrix(reduction),
      stdev = rep(0, ncol(reduction)),
      key = reduction.name,
      jackstraw = NULL,
      misc = list()))
    sc[[reduction.name]] <- DR
  } else if (inherits(sc, "SingleCellExperiment")) {
    reducedDim(sc, reduction.name) <- reduction
  }
  return(sc)
}

AF.col <- c(2,3,4,5,6)
KF.col <- c(7,8,9,10,11,12,13,14,15,16)


#Working with only barcodes that have CTaa
filter.object <- function(sc) {
  meta <- grabMeta(sc)
  cells.chains <- rownames(meta[!is.na(meta[["CTaa"]]),])
  if (inherits(x=sc, what ="Seurat")) {
    sc <- subset(sc, cells = cells.chains)
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
    sc <- sc[,colnames(sc) %in% cells.chains]
  }
  return(sc)
}

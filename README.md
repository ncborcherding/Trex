# Trex
Using TCR and expression for graph embedding

<img align="right" src="https://github.com/ncborcherding/Trex/blob/main/www/trex_hex.png" width="305" height="352">

## Introduction
Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/ncborcherding/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. However, utilization of clonal indices for more complex analyses are are still lacking, spefically in using clonality in embedding of single-cells and in trajectory analyses. To this end, here we develop the basis of combining clonal and expression analyses to facilitate dimensional reduction and trajectory inference. 

If you are looking for the (very cool) TCR-epitope prediction algorithm **TCRex**, check out their website [here](https://tcrex.biodatamining.be/).

# System requirements 

Trex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have had TCRs added using [scRepertoire](https://github.com/ncborcherding/scRepertoire). Trex has been tested on OS X and Windows platforms.

**muxViz** is essiential before installing Trex:

```r
devtools::install_github("manlius/muxViz")
```

**keras** is nessecary to use the autoencoder function (this includes the set up of the tensorflow environment in R):

```r
##Install keras
install.packages("keras")

##Setting up Tensor Flow
library(reticulate)
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(tensorflow)
install_tensorflow()
```

# Installation

To run Trex, open R and install Trex from github: 

```r
devtools::install_github("ncborcherding/Trex")
```

# Usage/Demos

Trex should be able to be run in popular R-based single-cell workflows, including Seurat and Bioconductor/Single-Cell Experiment formats.

## Quick Start 

Check out this [vignette](https://ncborcherding.github.io/vignettes/Trex.html) for a quick start tutorial. 

## Eigen Vector Matrix

The Trex algorithm allows users to select TCR-based metrics to return eigen vector to be used in dimensional reduction. If single-cell objects are not filtered for T cells with TCR,  `maTrex()` will still return values, however TREX_1 will be based on the disparity of TCR-containing and TCR-non-containg cells based on the Trex algorithm. 

```r
library(Trex)
my_trex <- maTrex(singleObject)
```

## Seurat or Single-Cell Experiment

You can run Trex within your Seurat or SingleCellExperiemt workflow. **Importantly** `runTrex()` will automatically filter single-cells that do not contain TCR information in the meta data of the single-cell object. 

```r
seuratObj_Tonly <- runTrex(seuratObj, #The single cell object
                   chains = "both", #Use of "TRA", "TRB" or "both"
                   edit.method = "lv", #Calculate edit distance methods
                   AA.properties = c("AF", "KF", "both"), 
                   AA.method = "auto", #Use "auto" for Autoencoder or 
                   #"mean" for mean properties across cdr3 sequence
                   reduction.name = "Trex", #Name designation for 
                   #the vectors to be added to the single-cell object
                   c.trim = 0, #Amino acid residues to trim from the start of the cdr3 sequence
                   n.trim = 0, #Amino acid residues to trim from the end of the cdr3 sequence
                   nearest.method = "threshold", #Use "threshold" to find related clonotypes above threshold or 
                   #"nn" for nearest neighbor by normalized distances
                   threshold = 0.85, #The normalized threshold to use where 1 = clone
                   near.neighbor = NULL, #The number of nearest neighbors to 
                   #use in generating an adjacency matrix
                   add.INKT = TRUE, #Add a additional layer for invariant natural 
                   #killer T cells based on genes
                   add.MAIT = TRUE, #Add a additional layer for Mucosal-associated 
                   #invariant T cells based on genes
                   species = "human") #Indicate "human" or "mouse" for gene-based metrics
                   
seuratObj_Tonly <- runTrex(seuratObj, reduction.name = "Trex")
```

## After Running Trex

From here, you can generate a tSNE/UMAP using the Trex values, similar to the PCA values based on variable gene expression.

```r
seuratObj <- RunTSNE(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
seuratObj <- RunUMAP(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
```

If using Seurat package, the Trex embedding infromation and gene expression PCA can be used to find the [Weighted Nearest Neighbors](https://pubmed.ncbi.nlm.nih.gov/34062119/). Before applying the WNN approach, best practice would be to remove the TCR-related genes from the list of variable genes and rerunning the PCA analysis. 

### Recaluclate PCA without TCR genes with queitTCRgenes() function in Trex.
```r
seuratObj <- quietTCRgenes(seuratObj)
seuratObj <- RunPCA(seuratObj)
```

### Running WNN approach
```r
seuratObj <- FindMultiModalNeighbors(
  seuratObj, reduction.list = list("pca", "Trex"), 
  dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight"
)
seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
```
***

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

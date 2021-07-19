# Trex
Using TCR and expression for graph embedding

<img align="right" src="https://github.com/ncborcherding/Trex/blob/main/www/trex_hex.png" width="305" height="352">

### Introduction
Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/ncborcherding/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. However, utilization of clonal indices for more complex analyses are are still lacking, spefically in using clonality in embedding of single-cells and in trajectory analyses. To this end, here we develop the basis of combining clonal and expression analyses to facilitate dimensional reduction and trajectory inference. 

# System requirements 

Trex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have had TCRs added using[scRepertoire](https://github.com/ncborcherding/scRepertoire). Trex has been tested on OS X and Windows platforms.

# Installation

To run Trex, open R and install harmony from CRAN: 

```r
library(devtools)
install_github("ncborcherding/Trex")
```

# Usage/Demos

Trex should be able to be run in popular R-based single-cell workflows, including Seurat and Bioconductor formats.

## Quick Start 

Check out this [vignette](https://ncborcherding.github.io/vignettes/Trex.html) for a quick start tutorial. 

## Eigen value matrix

The Trex algorithm iteratively corrects PCA embeddings. To input your own low dimensional embeddings directly, set `do_pca=FALSE`. Harmony is packaged with a small dataset 

```r
library(Trex)
my_trex <- calculatemaTrex(singleObject)
```

## Seurat 

You can run Trex within your Seurat or SingleCellExperiemt workflow. 

```r
seuratObj <- runTrex(seuratObj, "Trex")
```

From here, you can generate a tSNE/UMAP using the Trex values, similar to the PCA values based on variable gene expression.

```r
seuratObj <- RunTSNE(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
seuratObj <- RunUMAP(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
```

If using Seurat package, the Trex embedding infromation and gene expression PCA can be used to find the [Weighted Nearest Neighbors](https://pubmed.ncbi.nlm.nih.gov/34062119/)

```r
seuratObj <- FindMultiModalNeighbors(
  seuratObj, reduction.list = list("pca", "Trex"), 
  dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight"
)
seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
```

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

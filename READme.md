# GraphST-R

**An R implementation of GraphST for spatial transcriptomics: identifying spatial domains from 10x Visium data with graph self-supervised contrastive learning.**

[GraphST](https://github.com/JinmiaoChenLab/GraphST) (Long et al., *Nature Communications*, 2023) is a widely used method for spatial clustering, data integration, and cell-type deconvolution. It only exists in Python, though, and much of the spatial transcriptomics community works in R with Seurat and Bioconductor. GraphST-R re-architects the method in R so those researchers can use it without leaving their existing workflow.

This project was developed as part of the I-ACED research program at Prairie View A&M University.

---

## What it does

- **Builds a spatial neighbor graph** from spot coordinates (k-nearest neighbors) and normalizes the adjacency matrix
- **Learns spot embeddings** with a graph convolutional encoder trained using self-supervised contrastive learning (a corrupted feature view plus a readout/discriminator)
- **Clusters spots into spatial domains** with `mclust`, Leiden, or Louvain, including an optional spatial refinement step
- **Maps scRNA-seq to spatial data (deconvolution)** by training a mapping matrix that projects cell types onto spots (in progress)

## Architecture

The core is an object-oriented `GraphST` class built with **R6**. Its neural network modules are written in **torch for R**, and **reticulate** handles interoperability with Python tools such as scanpy.

| File | Purpose |
| --- | --- |
| `GraphR.R` | `GraphST` R6 class: initialization, `train()`, `train_sc()`, `train_map()`, and loss functions (contrastive and noise cross-entropy) |
| `model.R` | torch modules: graph encoder, discriminator, readout (`AvgReadout`), and the scRNA-seq encoder |
| `Preprocess.R` | Gene filtering, HVG preprocessing, neighbor graph construction, feature permutation, contrastive labels, adjacency normalization, sparse-to-torch tensor conversion, seeding |
| `utils.R` | Clustering (`mclust`, Leiden, Louvain), resolution search, label refinement, and cell-to-spot projection |
| `Test1.Rmd` | End-to-end walkthrough notebook |
| `data/GraphData/` | Example 10x Visium dataset |

## Example dataset

The walkthrough uses sample **151673** from the human dorsolateral prefrontal cortex (DLPFC) 10x Visium dataset (Maynard et al., 2021). Its cortical layers and white matter are manually annotated, which makes it a standard benchmark for spatial domain identification (7 clusters).

## Getting started

### Requirements

- R (≥ 4.2) and RStudio
- R packages:

```r
install.packages(c("R6", "Matrix", "dplyr", "torch", "pbapply", "Dict", "transport",
                   "Seurat", "SeuratObject", "FNN", "reticulate", "igraph",
                   "gmodels", "mclust", "progress"))
torch::install_torch()

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("SpatialExperiment", "DESeq2"))
```

- A Python environment with `scanpy`, used through reticulate

### Run the walkthrough

1. Clone the repo and open `ZoranayFirstProject.Rproj` in RStudio.
2. Unzip `data/GraphData.zip` if `data/GraphData/outs` isn't already present.
3. In `Test1.Rmd`, point `use_condaenv()` at your own Python environment.
4. Run the notebook:

```r
source("Preprocess.R")
source("utils.R")
source("model.R")
source("GraphR.R")

spe <- Load10X_Spatial(data.dir = "data/GraphData/outs")

model   <- GraphST$new(adata = spe)
trained <- model$train()

trained <- clustering(trained, n_clusters = 7, radius = 50,
                      method = "mclust", refinement = TRUE)
```

## Project status

This is a **research prototype**. Spatial domain identification is the main working path. Visualization and the scRNA-seq deconvolution workflow are still in progress, and the code is not yet packaged as an installable R package.

**Roadmap:**
- [ ] Package as an installable R package (`DESCRIPTION`, `NAMESPACE`, roxygen docs)
- [ ] Add spatial domain plots and ARI scoring against the manual DLPFC annotations
- [ ] Finish the deconvolution workflow
- [ ] Add unit tests (`testthat`)

## Acknowledgments

- Original method: Long, Y. et al. *Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST.* Nature Communications 14, 1155 (2023). [Python implementation](https://github.com/JinmiaoChenLab/GraphST)
- Dataset: Maynard, K.R. et al. *Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.* Nature Neuroscience (2021)
- I-ACED Research Program, Prairie View A&M University

## Author

**Zoranay Coursey** · [LinkedIn](https://www.linkedin.com/in/zoranay-coursey/) · [GitHub](https://github.com/myanna1416)

## License

Released under the [MIT License](LICENSE).





# load libraries 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)
# library(babelgene)

# Load processed seurat object
seu_rpca <- readRDS("03-analysis_scratch/seu_v5_preintegration.rds")

# Integrate with RPCA
seu_rpca <- IntegrateLayers(
  object = seu_rpca,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  assay = "RNA",
  verbose = FALSE
)

seu_rpca <- FindNeighbors(
  seu_rpca,
  reduction = "integrated.rpca",
  dims = 1:20,
  graph.name = "rpca_snn"
)

seu_rpca <- FindClusters(
  seu_rpca,
  graph.name = "rpca_snn",
  resolution = c(0.2,0.3,0.4)
)

seu_rpca <- RunUMAP(
  seu_rpca,
  reduction = "integrated.rpca",
  dims = 1:20,
  reduction.name = "umap.rpca"
)



DimPlot(
  seu_rpca,
  reduction = "umap.rpca",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

DimPlot(
  seu_rpca,
  reduction = "umap.rpca",
  group.by = "rpca_snn_res.0.3",
  split.by = "condition",
  label = TRUE
)


# use RNA counts/normalized data for DE
DefaultAssay(seu_rpca) <- "RNA"
seu_rpca <- JoinLayers(seu_rpca, assay = "RNA")

Idents(seu_rpca) <- "rpca_snn_res.0.3"

saveRDS(seu_rpca, file = "03-analysis_scratch/seu_v5_rpca.rds")

# seu_rpca <- readRDS("03-analysis_scratch/seu_v5_rpca.rds")

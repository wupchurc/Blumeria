# load libraries 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)

library(presto)
library(scCustomize)
library(SummarizedExperiment)
library(RColorBrewer)
library(circlize)
library(tidyr)


# ----Create UMAPs with RPCA Integration----

seu_rpca <- readRDS(
  "03-analysis_scratch/seu_v5_preintegration.rds"
)

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
  cluster.name = "rpca_clusters",
  resolution = 0.2
)

seu_rpca <- RunUMAP(
  seu_rpca,
  reduction = "integrated.rpca",
  dims = 1:20,
  reduction.name = "umap.rpca"
)

saveRDS(
  seu_rpca,
  file = "03-analysis_scratch/seu_v5_rpca.rds"
)

seu_rpca <- readRDS("03-analysis_scratch/seu_v5_rpca.rds")

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
  group.by = "rpca_clusters",
  split.by = "condition",
  label = TRUE
)

#----Create UMAPs with Harmony Integration----

seu_harmony <- readRDS(
  "03-analysis_scratch/seu_v5_preintegration.rds"
)

seu_harmony <- IntegrateLayers(
  object = seu_harmony,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  assay = "RNA",
  verbose = FALSE
)

seu_harmony <- FindNeighbors(
  seu_harmony,
  reduction = "harmony",
  dims = 1:20,
  graph.name = "harmony_snn"
)

seu_harmony <- FindClusters(
  seu_harmony,
  graph.name = "harmony_snn",
  cluster.name = "harmony_clusters",
  resolution = 0.3
)

seu_harmony <- RunUMAP(
  seu_harmony,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap.harmony"
)

saveRDS(
  seu_harmony,
  file = "03-analysis_scratch/seu_v5_harmony.rds"
)

seu_harmony <- readRDS("03-analysis_scratch/seu_v5_harmony.rds")


DimPlot(
  seu_harmony,
  reduction = "umap.harmony",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

DimPlot(
  seu_harmony,
  reduction = "umap.harmony",
  group.by = "harmony_clusters",
  split.by = "condition",
  label = TRUE
)
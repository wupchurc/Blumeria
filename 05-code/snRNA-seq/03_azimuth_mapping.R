library(Azimuth)

# 1. Convert rat Seurat object features to human orthologs
seu_human <- orthogene::convert_orthologs(
  gene_df = seu_rpca,
  input_species = "rat",
  output_species = "human"
)

# 2. Run Azimuth against the Human Heart Atlas reference
seu_human <- RunAzimuth(seu_human, reference = "heartref")

# 3. Transfer predicted annotations back to your original rat object
seu_rpca$azimuth_celltype <- seu_human$predicted.annotation.l2

# View automated predictions on UMAP
DimPlot(seu_rpca, reduction = "umap.rpca", group.by = "azimuth_celltype", label = TRUE)
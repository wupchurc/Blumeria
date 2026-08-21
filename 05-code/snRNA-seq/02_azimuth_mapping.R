# ==============================================================================
# Script: 02_azimuth_mapping.R
# Purpose: Convert rat Seurat object to human orthologs and run Azimuth 
#          mapping against the Human Heart Atlas reference to automate
#          cell-type predictions.
# Input:   03-analysis_scratch/seu_v5_rpca.rds
# Output:  03-analysis_scratch/seu_v5_azimuth.rds
# ==============================================================================

# 1. Load Required Libraries ----
library(Seurat)
library(Azimuth)
library(orthogene)
library(patchwork)
library(ggplot2)

# 2. Load Integrated Seurat Object ----
message("Loading integrated Seurat object...")
seu_rpca <- readRDS("03-analysis_scratch/seu_v5_rpca.rds")

# Ensure default assay is set to RNA for proper feature extraction
DefaultAssay(seu_rpca) <- "RNA"

# 3. Convert Species: Rat to Human ----
message("Extracting raw rat count matrix...")
# Extract the raw count data rather than passing the whole Seurat object
rat_counts <- LayerData(seu_rpca, assay = "RNA", layer = "counts")

message("Converting rat gene features to human orthologs...")
# Convert the matrix row names from rat to human
human_counts <- orthogene::convert_orthologs(
  gene_df = rat_counts,
  gene_input = "rownames", 
  input_species = "rat",
  output_species = "human",
  method = "gprofiler",
  non121_strategy = "drop_both_species" # Safely drops ambiguous multi-mapped genes
)

message("Creating temporary human Seurat object...")
# Wrap the converted matrix back into a barebones Seurat object for Azimuth
seu_human <- CreateSeuratObject(counts = human_counts)

# Verify conversion (ensure row names look like human HGNC symbols, e.g., "TRIM63")
print(head(rownames(seu_human)))

# 4. Run Azimuth Reference Mapping ----
# We use "heartref" corresponding to the Human Heart Cell Atlas reference
message("Running Azimuth mapping against 'heartref'...")
seu_human <- RunAzimuth(
  query = seu_human,
  reference = "heartref"
)

# 5. Transfer Annotations Back to Original Object ----
# The heartref uses 'celltype' instead of 'annotation' for its metadata columns
message("Transferring predictions back to the rat Seurat object...")
seu_rpca$azimuth_l1 <- seu_human$predicted.celltype.l1
seu_rpca$azimuth_l2 <- seu_human$predicted.celltype.l2
seu_rpca$azimuth_score <- seu_human$predicted.celltype.l2.score

# 6. Visualize Predictions ----
# Plot original unsupervised clusters vs Azimuth predictions side-by-side
p1 <- DimPlot(seu_rpca, reduction = "umap.rpca", group.by = "rpca_snn_res.0.3", label = TRUE) + 
  ggtitle("Unsupervised Clusters (Res 0.3)")

p2 <- DimPlot(seu_rpca, reduction = "umap.rpca", group.by = "azimuth_l2", label = TRUE, repel = TRUE) + 
  ggtitle("Azimuth Predicted Cell Types (Level 2)")

# Display plot
p1 | p2

# Optional: Save visualization to file
ggsave("03-analysis_scratch/umap_azimuth_predictions.pdf", plot = (p1 | p2), width = 20, height = 6)

# 7. Save the Updated Object ----
message("Saving updated Seurat object...")
saveRDS(seu_rpca, file = "03-analysis_scratch/seu_v5_azimuth.rds")
message("Done.")
# ==============================================================================
# Script: 03_markers_annotation.R
# Purpose: Manually query individual unsupervised clusters to extract top 
#          rat marker genes for biological verification of Azimuth labels.
# Input:   03-analysis_scratch/seu_v5_azimuth.rds
# Output:  Targeted marker lists and (eventually) annotated Seurat object
# ==============================================================================

# 1. Load Required Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. Load the Azimuth-Mapped Object ----
message("Loading Azimuth-mapped Seurat object...")
seu_rpca <- readRDS("03-analysis_scratch/seu_v5_azimuth.rds")

# Ensure default assay is RNA for accurate raw expression testing
DefaultAssay(seu_rpca) <- "RNA"

# Set the active identity to your unsupervised clusters
Idents(seu_rpca) <- "rpca_snn_res.0.6"

# 3. Manually Query a Specific Cluster for Conserved Markers ----
target_cluster <- "0"

message(paste("Finding conserved markers for Cluster", target_cluster, "..."))
cluster_markers <- FindConservedMarkers(
  object = seu_rpca,
  ident.1 = target_cluster,
  grouping.var = "condition", # Ensures markers are present in all experimental groups
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# 4. View and Save the Top Markers ----
# FindConservedMarkers automatically calculates a 'max_pval' across your groups.
# We filter for genes that are significantly upregulated in ALL conditions.
top_markers <- cluster_markers %>%
  filter(max_pval < 0.05) 

# Print the top 20 genes to the console. 
# Note: The output will have columns for each of your conditions (e.g., Control_avg_log2FC)
print(head(top_markers, 100))
cat(paste(rownames(head(top_markers, 20)), collapse = ", "), "\n")

# ==============================================================================
# 5. Master Cluster Annotation (Fill this out sequentially!) ----
# ==============================================================================
# Reset identities to the original numbers just to be safe before renaming
Idents(seu_rpca) <- "rpca_snn_res.0.6"

# Build your dictionary as you identify each cluster
seu_rpca <- RenameIdents(
  seu_rpca,
  "0" = "Fibroblast",
  "1" = "Fibroblast",
  "2" = "Pericyte",
  "3" = "Vascular EC",
  "4" = "Vascular EC",
  "5" = "Cardiomyocyte",
  "6" = "Cardiomyocyte",
  "7" = "Macrophage",
  "8" = "Vascular EC",
  "9" = "Macrophage",
  "10" = "Fibroblast",
  "11" = "Cardiomyocyte",
  "12" = "T Cell",
  "13" = "Endocardial EC",
  "14" = "NK Cell",
  "15" = "Fibroblast",
  "16" = "Dendritic Cell",
  "17" = "Monocyte",
  "18" = "Lymphatic EC",
  "19" = "B Cell",
  "20" = "Neuronal", # Glial Cell
  "21" = "Neutrophil",
  "22" = "Vascular EC",
  "23" = "Epicardial Cell",
  "24" = "Doublet",
  "25" = "T Cell",
  "26" = "Pericyte" # or smooth muscle cell
)

# 6. Save Annotations to Metadata and Export ----
# Storing this in a new metadata column ensures you never lose your original numbers
seu_rpca$manual_annotation <- Idents(seu_rpca)

# Subset to remove doublets for clean downstream comparison
seu_clean <- subset(seu_rpca, idents = "Doublet", invert = TRUE)

# Save your fully annotated object
saveRDS(seu_clean, "03-analysis_scratch/seu_annotated.rds")

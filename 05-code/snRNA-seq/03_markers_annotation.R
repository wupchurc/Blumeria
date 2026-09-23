# ==============================================================================
# Script: 03_markers_annotation.R
# Purpose: Manually query individual unsupervised clusters to extract top 
#          rat marker genes for biological verification of Azimuth labels.
# Input:   03-analysis_scratch/seu_v5_azimuth.rds
# Output:  Targeted marker lists and annotated Seurat object
#          Heatmap of cell markers
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
top_markers <- cluster_markers %>% filter(max_pval < 0.05) 

# Print the top genes to the console. 
# Note: The output will have columns for each of your conditions (e.g., Control_avg_log2FC)
print(head(top_markers, 100))
cat(paste(rownames(head(top_markers, 20)), collapse = ", "), "\n")

# ==============================================================================
# Master Cluster Annotation ----
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
  "23" = "Doublet",                #"Epicardial Cell",
  "24" = "Doublet",
  "25" = "T Cell",
  "26" = "Pericyte" # or smooth muscle cell
)

# ==============================================================================
# Save Annotations & Order Metadata ----
# ==============================================================================
# Save annotations to metadata slot
seu_rpca$manual_annotation <- Idents(seu_rpca)

# Subset to remove doublets for clean downstream comparison
seu_clean <- subset(seu_rpca, idents = "Doublet", invert = TRUE)

# Set clean factor levels
# celltype_order <- c(
  # "Cardiomyocyte", "Fibroblast", "Pericyte", "Vascular EC", 
  # "Endocardial EC", "Lymphatic EC", "Macrophage", 
  # "Monocyte", "Dendritic Cell", "Neutrophil", 
  # "T Cell", "B Cell", "NK Cell", "Neuronal"
# )

celltype_order <- c(
  "Cardiomyocyte", "Fibroblast", "Macrophage", "Monocyte",
  "Dendritic Cell", "Neutrophil", "T Cell", "B Cell", "NK Cell",
  "Pericyte", "Vascular EC", "Endocardial EC", "Lymphatic EC",
  "Neuronal"
)

# Drop unused levels and order factor
existing_levels <- intersect(celltype_order, unique(seu_clean$manual_annotation))
seu_clean$manual_annotation <- factor(seu_clean$manual_annotation, levels = existing_levels)
Idents(seu_clean) <- "manual_annotation"

# Save cleaned and ordered object
saveRDS(seu_clean, "03-analysis_scratch/seu_annotated.rds")

# ==============================================================================
# 7. Generate Marker Heatmap for Annotation Validation ----
# ==============================================================================

# 1. Define Aesthetics (Matching Script 04 Palette and Order) ----
cell_cols <- c(
  "Cardiomyocyte"  = "#BA1C30", "Fibroblast"     = "#5FA641", 
  "Macrophage"     = "#702C8C", "Monocyte"       = "#CC79A7",
  "Dendritic Cell" = "#D06B70", "Neutrophil"     = "#D55E00",
  "T Cell"         = "#999999", "B Cell"         = "#4DAF4A", 
  "NK Cell"        = "#377EB8", "Pericyte"       = "#0072B2",
  "Vascular EC"    = "#F0E442", "Endocardial EC" = "#009E73",
  "Lymphatic EC"   = "#56B4E9", "Neuronal"       = "#E69F00"
)

cell_abbr <- c(
  "Cardiomyocyte"  = "CM",      "Fibroblast"     = "FB",
  "Macrophage"     = "M\u03A6", "Monocyte"       = "Mono",
  "Dendritic Cell" = "DC",      "Neutrophil"     = "Neu",
  "T Cell"         = "T",       "B Cell"         = "B",
  "NK Cell"        = "NK",      "Pericyte"       = "Peri",
  "Vascular EC"    = "V EC",     "Endocardial EC" = "E EC",
  "Lymphatic EC"   = "L EC",     "Neuronal"       = "Neur"
)

# Align colors with abbreviated names
cell_cols_abbr <- setNames(cell_cols, cell_abbr[names(cell_cols)])

# 2. Add Abbreviated Metadata & Set Cluster Order ----
seu_clean$cell_type_abbr <- unname(cell_abbr[as.character(seu_clean$manual_annotation)])
abbr_order <- unname(cell_abbr[celltype_order])

# Retain only levels existing in dataset
existing_abbr_levels <- intersect(abbr_order, unique(seu_clean$cell_type_abbr))
seu_clean$cell_type_abbr <- factor(seu_clean$cell_type_abbr, levels = existing_abbr_levels)
Idents(seu_clean) <- "cell_type_abbr"

# 3. Define & Scale Canonical Markers ----
DefaultAssay(seu_clean) <- "RNA"

validation_markers <- c(
  "Tnnt2", "Myh6", "Ryr2",     # Cardiomyocyte
  "Pdgfra", "Col1a1", "Dcn",   # Fibroblast
  "Mrc1", "Mertk", "C1qa",     # Macrophage
  "Plac8", "Itgal", "Spn",     # Monocyte
  "Flt3", "Ciita", "Wdfy4",    # Dendritic Cell
  "S100a8", "S100a9", "Csf3r", # Neutrophil
  "Cd3e", "Cd4", "Cd8a",       # T Cell
  "Ms4a1", "Cd79a", "Pax5",    # B Cell
  "Ncr1", "Klrk1",             # NK Cell
  "Kcnj8", "Pdgfrb", "Rgs5",   # Pericyte
  "Fabp4", "Aqp1", "Sox17",    # Vascular EC
  "Npr3", "Pkhd1l1", "Nrg1",   # Endocardial EC
  "Prox1", "Pdpn", "Ccl21",    # Lymphatic EC
  "Cdh19", "Scn7a", "Lgi4"     # Neuronal
)

valid_features <- intersect(validation_markers, rownames(seu_clean))
seu_clean <- ScaleData(seu_clean, features = valid_features, verbose = FALSE)

# 4. Downsample & Generate Heatmap ----
seu_heatmap_sub <- subset(seu_clean, downsample = 300)

p_heatmap <- DoHeatmap(
  object           = seu_heatmap_sub,
  features         = valid_features,
  group.by         = "cell_type_abbr",
  group.colors     = cell_cols_abbr,
  size             = 4,
  angle            = 0,            # Keeps text horizontal inside top bar
  hjust            = 0.5,          # Centers text horizontally
  vjust            = -1.75,          # Centers text vertically in the color bar
  group.bar.height = 0.05,
  draw.lines       = TRUE,
  disp.min         = -2.5, 
  disp.max         = 2.5
) + 
  guides(color = "none") +  
  scale_fill_gradient2(
    low      = "magenta", 
    mid      = "black", 
    high     = "yellow", 
    midpoint = 0, 
    limits   = c(-2.5, 2.5),
    name     = "Z-score"
  ) +
  theme(
    axis.text.y = element_text(size = 9, face = "italic")
  )

print(p_heatmap)

# Save figure
ggsave(
  filename = "04-results/annotation_validation_heatmap.png",
  plot     = p_heatmap,
  width    = 12,
  height   = 9,
  units    = "in",
  dpi      = 300
)

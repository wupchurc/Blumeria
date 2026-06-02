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

# load processed seurat object from rds file 
seu_integrated <- readRDS("03-analysis_scratch/seu_integrated_strict.rds")

seu_integrated <- FindNeighbors(seu_integrated, reduction = "harmony", dims = 1:20)
seu_integrated <- FindClusters(seu_integrated, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
seu_integrated <- RunUMAP(seu_integrated, reduction = "harmony", dims = 1:20)

seu_integrated$condition <- factor(
  seu_integrated$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)
DimPlot(
  seu_integrated,
  reduction = "umap",
  group.by = "condition",
  split.by = "condition",
  label = FALSE
)

DimPlot(
  seu_integrated,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  split.by = "condition",
  label = TRUE
)

# check percent.mt per cluster
VlnPlot(seu_integrated, features = c("percent.mt"), group.by = "RNA_snn_res.0.3") 
  # geom_hline(yintercept = 15, linetype = "dashed", color = "red")
  

# use RNA counts/normalized data for DE
DefaultAssay(seu_integrated) <- "RNA"
seu_integrated <- JoinLayers(seu_integrated, assay = "RNA")

Idents(seu_integrated) <- "RNA_snn_res.0.3"

# ---- conserved markers per cluster across conditions ----
# Run FindConservedMarkers on clusters to help identify
markers_cluster0  <- FindConservedMarkers(seu_integrated, ident.1 = 0,  grouping.var = "condition")
markers_cluster14 <- FindConservedMarkers(seu_integrated, ident.1 = 14, grouping.var = "condition")
markers_cluster4 <- FindConservedMarkers(seu_integrated, ident.1 = 4, grouping.var = "condition")
markers_cluster8 <- FindConservedMarkers(seu_integrated, ident.1 = 8, grouping.var = "condition")
markers_cluster12 <- FindConservedMarkers(seu_integrated, ident.1 = 12, grouping.var = "condition")
markers_cluster13 <- FindConservedMarkers(seu_integrated, ident.1 = 13, grouping.var = "condition")
markers_cluster2 <- FindConservedMarkers(seu_integrated, ident.1 = 2, grouping.var = "condition")
markers_cluster5 <- FindConservedMarkers(seu_integrated, ident.1 = 5, grouping.var = "condition")
markers_cluster26 <- FindConservedMarkers(seu_integrated, ident.1 = 26, grouping.var = "condition")
markers_cluster1 <- FindConservedMarkers(seu_integrated, ident.1 = 1, grouping.var = "condition")
markers_cluster6 <- FindConservedMarkers(seu_integrated, ident.1 = 6, grouping.var = "condition")
markers_cluster3 <- FindConservedMarkers(seu_integrated, ident.1 = 3, grouping.var = "condition")
markers_cluster9 <- FindConservedMarkers(seu_integrated, ident.1 = 9, grouping.var = "condition")
markers_cluster11 <- FindConservedMarkers(seu_integrated, ident.1 = 11, grouping.var = "condition")
markers_cluster10 <- FindConservedMarkers(seu_integrated, ident.1 = 10, grouping.var = "condition")
markers_cluster7 <- FindConservedMarkers(seu_integrated, ident.1 = 7, grouping.var = "condition")
markers_cluster16 <- FindConservedMarkers(seu_integrated, ident.1 = 16, grouping.var = "condition")
markers_cluster18 <- FindConservedMarkers(seu_integrated, ident.1 = 18, grouping.var = "condition")
markers_cluster17 <- FindConservedMarkers(seu_integrated, ident.1 = 17, grouping.var = "condition")
markers_cluster22 <- FindConservedMarkers(seu_integrated, ident.1 = 22, grouping.var = "condition")
markers_cluster21 <- FindConservedMarkers(seu_integrated, ident.1 = 21, grouping.var = "condition")
markers_cluster23 <- FindConservedMarkers(seu_integrated, ident.1 = 23, grouping.var = "condition")
markers_cluster24 <- FindConservedMarkers(seu_integrated, ident.1 = 24, grouping.var = "condition")
markers_cluster25 <- FindConservedMarkers(seu_integrated, ident.1 = 25, grouping.var = "condition")
markers_cluster15 <- FindConservedMarkers(seu_integrated, ident.1 = 15, grouping.var = "condition")
markers_cluster19 <- FindConservedMarkers(seu_integrated, ident.1 = 19, grouping.var = "condition")
markers_cluster20 <- FindConservedMarkers(seu_integrated, ident.1 = 20, grouping.var = "condition")

# Renaming clusters 
print(Idents(seu_integrated))
seu_integrated <- RenameIdents(seu_integrated, '0' = 'Fibroblasts')
seu_integrated <- RenameIdents(seu_integrated, '14' = 'Neutrophils')
seu_integrated <- RenameIdents(seu_integrated, '4' = 'Pericytes')
seu_integrated <- RenameIdents(seu_integrated, '8' = 'Endocardial EC')
seu_integrated <- RenameIdents(seu_integrated, '12' = 'Lymphatic EC')
seu_integrated <- RenameIdents(seu_integrated, '13' = 'Lymphocytes')
seu_integrated <- RenameIdents(seu_integrated, '26' = 'Lymphocytes')
seu_integrated <- RenameIdents(seu_integrated, '2' = 'Vascular EC')
seu_integrated <- RenameIdents(seu_integrated, '5' = 'Vascular EC')
seu_integrated <- RenameIdents(seu_integrated, '1' = "Cardiomyocytes")
seu_integrated <- RenameIdents(seu_integrated, '6' = "Cardiomyocytes")
seu_integrated <- RenameIdents(seu_integrated, '3' = 'Macrophages')
seu_integrated <- RenameIdents(seu_integrated, '10' = 'Lymphocytes')
seu_integrated <- RenameIdents(seu_integrated, '7' = 'Lymphocytes')
seu_integrated <- RenameIdents(seu_integrated, '18' = 'Fibroblasts')
seu_integrated <- RenameIdents(seu_integrated, '22' = 'Pericytes')
seu_integrated <- RenameIdents(seu_integrated, '21' = 'Pericytes')
seu_integrated <- RenameIdents(seu_integrated, '23' = "Cardiomyocytes")
seu_integrated <- RenameIdents(seu_integrated, '15' = 'Glial Cells')
seu_integrated <- RenameIdents(seu_integrated, '9' = 'Monocytes')
seu_integrated <- RenameIdents(seu_integrated, '16' = 'Noise') # high Mt-co2, Mt-co3, Mt-cyb, Mt-nd4
seu_integrated <- RenameIdents(seu_integrated, '17' = 'Lymphocytes')
seu_integrated <- RenameIdents(seu_integrated, '11' = 'Dendritic Cells')
seu_integrated <- RenameIdents(seu_integrated, '24' = 'Noise') # contaminated immune cells?
seu_integrated <- RenameIdents(seu_integrated, '25' = 'Macrophages')
seu_integrated <- RenameIdents(seu_integrated, '19' = 'Monocytes')
seu_integrated <- RenameIdents(seu_integrated, '20' = 'Lymphocytes')


# Remove noise clusters
# List all clusters first
table(Idents(seu_integrated))

# Define clusters to REMOVE
bad_clusters <- c("Noise")

# Get all cells from bad clusters
bad_cells <- unlist(lapply(bad_clusters, function(x) WhichCells(seu_integrated, idents = x)))

# 4. Remove them ALL at once
seu_integrated <- seu_integrated[, !colnames(seu_integrated) %in% bad_cells]

# 5. Clean up clusters & refresh UMAP
# Idents(seu_integrated) <- "seurat_clusters"
# seu_integrated <- RunUMAP(seu_integrated, dims = 1:20)
DimPlot(seu_integrated, label = TRUE)
table(Idents(seu_integrated))  # Verify bad clusters gone 





# Store renamed Idents in metadata column
seu_integrated$cell_type <- Idents(seu_integrated)

# Set the desired order of cell types
celltype_order <- c(
  "Cardiomyocytes",
  "Fibroblasts",
  "Macrophages",
  "Monocytes",
  "Dendritic Cells",
  "Neutrophils",
  "Lymphocytes",
  "Pericytes",
  "Vascular EC",
  "Endocardial EC",
  "Lymphatic EC",
  # "Glial Cells"
  "Neural"
)

# Define colors matching celltype_order 
cell_cols <- c(
  "Cardiomyocytes" = "#BA1C30",
  "Fibroblasts"    = "#5FA641", 
  "Macrophages"      = "#702C8C",
  "Monocytes"      = "#CC79A7",
  "Dendritic Cells" = "#D06B70",
  "Neutrophils"   = "#D55E00",
  "Lymphocytes"   = "#999999",
  "Pericytes"    = "#0072B2",
  "Vascular EC"    = "#F0E442",
  "Endocardial EC"        = "#009E73",
  "Lymphatic EC"        = "#56B4E9",
  # "Glial Cells"       = "#E69F00"
  "Neural" = "#E69F00"
)

# Make cell_type a factor with that order
seu_integrated$cell_type <- factor(
  seu_integrated$cell_type,
  levels = celltype_order
)

# ---- Heatmap of markers used to validate cell type identities ----
all_markers <- c(
  "Tnnt2","Myh6","Ryr2", # cardiomyocyte markers
  "Pdgfra", "Col1a1", "Dcn", # fibroblast markers
  "Mrc1","Mertk","C1qa", #macrophage markers
  "Plac8","Itgal","Spn", # monocytes markers
  "Flt3","Ciita", "Wdfy4",            # dendritic cells
  "S100a8", "S100a9","Csf3r",              # neutrophil markers
  "Cd3e","Ms4a1","Ncr1",     # lymphocyte markers
  "Kcnj8","Pdgfrb","Rgs5", # pericyte markers
  "Fabp4","Aqp1","Sox17",# vascular ec markers
  "Npr3","Pkhd1l1","Nrg1",# endocardial ec markers
  "Prox1", "Pdpn", "Ccl21",   # lymphatic ec markers
  "Cdh19", "Scn7a", "Lgi4"                     # neuronal/ glia marker
)

p_heat <- DoHeatmap(subset(seu_integrated, downsample = 1000),
          features = all_markers,
          size = 3,
          hjust = 0.5,
          vjust = - 1.5,
          angle = 0, 
          group.bar.height = 0.06,
          group.by = "cell_type",
          group.colors = cell_cols) + 
  guides(color = "none") + 
  labs(fill = "Z-score") 

print(p_heat)

ggsave(
  filename = "04-results/celltype_markers_heatmap.png",  # or .png, .tiff
  plot     = p_heat,
  width    = 16,   # adjust as needed
  height   = 6,
  dpi      = 300
)

# --- plotting (still using integrated UMAP) ----
DefaultAssay(seu_integrated) <- "integrated"  # for plotting/clustering context

# Reorder condition for plots
seu_integrated$condition <- factor(
  seu_integrated$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

# Plot UMAPs with identified clusters
print(DimPlot(seu_integrated, reduction = 'umap', 
              label = TRUE, label.size = 4, repel = TRUE, cols = cell_cols) + NoLegend())

# ---- Relative quantification ----

# After setting celltype_order, reverse it for plotting consistency
plot_order <- rev(celltype_order)  # Reverse for ggplot/barplot
# Set Idents and cell_type with plot order
Idents(seu_integrated) <- factor(seu_integrated$cell_type, levels = plot_order)
seu_integrated$cell_type <- factor(seu_integrated$cell_type, levels = plot_order)
# Now cell_types and plots follow your desired order
cell_types <- levels(Idents(seu_integrated))

origin_x <- -15
origin_y <- -15
len      <- 5   # length of mini-axes in UMAP units

mini_axes <- data.frame(
  x    = c(origin_x, origin_x),
  y    = c(origin_y, origin_y),
  xend = c(origin_x + len, origin_x),
  yend = c(origin_y,        origin_y + len)
)
# Use same palette in UMAP
p_umap <- DimPlot(
  seu_integrated,
  reduction  = "umap",
  split.by   = "condition",
  label      = TRUE,
  # label      = FALSE,
  label.size = 4,
  cols       = cell_cols          # <— key line
) +
  NoGrid() + 
  NoLegend() +
  coord_cartesian(clip = "off") +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text   = element_blank(),
    axis.title  = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 30),
    legend.text = element_text(size = 20)
  ) +
  geom_segment(
    data = mini_axes,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(3, "pt")),
    linewidth = 0.4
  ) +
  annotate("text",
           x = origin_x + len, y = origin_y,
           label = "UMAP1", vjust = 1.5, size = 3) +
  annotate("text",
           x = origin_x - 0.25, y = origin_y + len,
           label = "UMAP2", hjust = 0.5, angle = 90, size = 3)

# counts per cell type per condition
tab_counts <- table(
  CellType  = Idents(seu_integrated),
  Condition = seu_integrated$condition
)

# relative proportions within each condition (columns sum to 1)
tab_prop <- prop.table(tab_counts, margin = 2)
df_prop <- as.data.frame(tab_prop) %>%
  dplyr::rename(Proportion = Freq)

# Extracting Percent Abundances
df_percentages_long <- df_prop %>%
  mutate(Percentage = round(Proportion * 100, 2)) %>%
  select(-Proportion) # Removes the fractional column
print("--- Long Format Percentages ---")
print(df_percentages_long)
df_percentages_wide <- df_percentages_long %>%
  pivot_wider(names_from = Condition, values_from = Percentage)
print("--- Wide Format Table (%) ---")
print(df_percentages_wide)
# Save to CSV
write.csv(df_percentages_wide, 
          file = "04-results/cell_type_percent_abundance.csv", 
          row.names = FALSE)

p_bar <- ggplot(df_prop,
                aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 1) +
  coord_flip() +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  xlab(NULL) +   # now appears above the bars
  ylab(NULL) +
  theme_minimal() +
  labs(title = "Relative Abundance (%)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),  # keep digits hidden
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),   
    axis.ticks.y     = element_blank()
  )

# plot combined graphs
p_combined <- p_umap / p_bar + plot_layout(heights = c(12, 1))

print(p_combined)





# ---- Create df_bar_data for labels ----

# 1. Define abbreviations to match your previous plot exactly
cell_abbr <- c(
  "Cardiomyocytes"  = "CM",
  "Fibroblasts"     = "FB",
  "Macrophages"     = "Mp", 
  "Monocytes"       = "Mono",
  "Dendritic Cells" = "DC",
  "Neutrophils"     = "Neu",
  "Lymphocytes"     = "Lym",
  "Pericytes"       = "Peri",
  "Vascular EC"     = "vEC",
  "Endocardial EC"  = "eEC",
  "Lymphatic EC"    = "lEC",
  "Neural"          = "Neural" # Mapping your new 'Neural' back to 'Glia' for the label
)

# 2. Calculate percentages, create the label string, and find midpoints
df_bar_data <- df_prop %>%
  mutate(
    # Keep the two decimal places you had previously
    Percentage = round(Proportion * 100, 2),
    
    # Combine Abbreviation, "=", Percentage, and "%" (e.g., "FB=35.28%")
    Label = paste0(cell_abbr[as.character(CellType)], "=", Percentage, "%") 
  ) %>%
  group_by(Condition) %>%
  # Reverse the CellType order to match ggplot's coord_flip stacking
  arrange(Condition, desc(CellType)) %>% 
  mutate(
    cumulative = cumsum(Proportion),
    midpoint = cumulative - (Proportion / 2)
  ) %>%
  ungroup()


# ---- 2. Include ALL, Evenly Space, and STAGGER labels tighter 
df_labels <- df_bar_data %>% 
  # Filter removed entirely to include all cell types
  group_by(Condition) %>%
  # Sort strictly left-to-right
  arrange(Condition, midpoint) %>% 
  mutate(
    # Spread the 12 labels evenly across the horizontal space
    # (Extended slightly to 0.88 to use the available width)
    label_pos = seq(0.0, 0.88, length.out = n()),
    
    # NEW: Bring the rows closer together
    # Top row at 0.55, Bottom row at 0.35 (less vertical gap)
    v_pos = ifelse(row_number() %% 2 != 0, 0.55, 0.35)
  ) %>%
  ungroup()


# ---- Build the modified Bar Plot 

p_bar <- ggplot(df_bar_data, aes(x = 1, y = Proportion, fill = CellType)) +
  # 1. The stacked bar
  geom_col(position = "stack", color = "black", width = 0.4) + 
  
  # 2. The colored circles (Reduced size to 2.5)
  geom_point(data = df_labels, aes(x = v_pos, y = label_pos, color = CellType), size = 2.5) +
  
  # 3. The text labels (Reduced size to 2.8, smaller nudge)
  geom_text(data = df_labels, aes(x = v_pos, y = label_pos, label = Label), 
            hjust = 0, nudge_y = 0.015, size = 2.8) +
  
  # Adjusted xlim to c(0.2, 1.2) to trim dead space at the bottom now that rows are tighter
  coord_flip(xlim = c(0.2, 1.2), ylim = c(0, 1), clip = "off") + 
  facet_wrap(~ Condition, nrow = 1) +
  
  scale_fill_manual(values = cell_cols, guide = "none") +
  scale_color_manual(values = cell_cols, guide = "none") +
  
  labs(title = "Relative Abundance (%)", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank()
  )

# plot combined graphs (Reduced height slightly from 2.5 to 2.2 since we tightened the rows)
p_combined <- p_umap / p_bar + plot_layout(heights = c(12, 2.2))

print(p_combined)

# ---- Save  ----

ggsave(
  filename = "04-results/umap_rel_abundance.png",
  plot     = p_combined,
  width    = 24,    # adjust as needed
  height   = 8,
  units = "in",
  dpi      = 300
)

saveRDS(seu_integrated, file = "03-analysis_scratch/seu_for_DGE.rds")

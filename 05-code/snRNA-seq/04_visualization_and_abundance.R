# ==============================================================================
# Script: 04_visualization_and_abundance.R
# Purpose: Generate publication-ready UMAPs and stacked bar plots of relative 
#          abundances across experimental conditions.
# Input:   03-analysis_scratch/seu_annotated.rds
# Output:  04-results/umap_rel_abundance.png
#          04-results/umap_combined.png
#          03-analysis_scratch/seu_after_visualization.rds
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# 1. Load Data & Ensure Ordering
seu_clean <- readRDS("03-analysis_scratch/seu_annotated.rds")

# Lock the condition levels so the facets and bars are ordered correctly
seu_clean$condition <- factor(
  seu_clean$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

# Initialize standard cell types from manual annotations
Idents(seu_clean) <- "manual_annotation"
seu_clean$cell_type <- Idents(seu_clean)

# 2. Define Plot Aesthetics (Order, Colors, Abbreviations)
celltype_order <- c(
  "Cardiomyocyte", "Fibroblast", "Macrophage", "Monocyte",
  "Dendritic Cell", "Neutrophil", "T Cell", "B Cell", "NK Cell",
  "Pericyte", "Vascular EC", "Endocardial EC", "Lymphatic EC",
  "Neuronal", "Epicardial Cell"
)

cell_cols <- c(
  "Cardiomyocyte"  = "#BA1C30", "Fibroblast"     = "#5FA641", 
  "Macrophage"     = "#702C8C", "Monocyte"       = "#CC79A7",
  "Dendritic Cell" = "#D06B70", "Neutrophil"     = "#D55E00",
  "T Cell"         = "#999999", "B Cell"         = "#4DAF4A", 
  "NK Cell"        = "#377EB8", "Pericyte"       = "#0072B2",
  "Vascular EC"    = "#F0E442", "Endocardial EC" = "#009E73",
  "Lymphatic EC"   = "#56B4E9", "Neuronal"       = "#E69F00",
  "Epicardial Cell"= "#A65628"  
)

cell_abbr <- c(
  "Cardiomyocyte"  = "CM",      "Fibroblast"     = "FB",
  "Macrophage"     = "M\u03A6", "Monocyte"       = "Mono",
  "Dendritic Cell" = "DC",      "Neutrophil"     = "Neu",
  "T Cell"         = "T",       "B Cell"         = "B",
  "NK Cell"        = "NK",      "Pericyte"       = "Peri",
  "Vascular EC"    = "vEC",     "Endocardial EC" = "eEC",
  "Lymphatic EC"   = "lEC",     "Neuronal"       = "Neur",
  "Epicardial Cell"= "Epi"
)

# Apply plotting order (reverse standard order so they stack correctly in ggplot)
plot_order <- rev(celltype_order)  

# 3. Clean Up Outlier Cells Prior to Metadata Assignments
# Convert the cell_type column to character to allow reassignment
seu_clean$cell_type <- as.character(seu_clean$cell_type)

# Locate the single stray Epicardial cell in MCT-Blumeria
stray_epi_idx <- seu_clean$cell_type == "Epicardial Cell" & seu_clean$condition == "MCT-Blumeria"

# Reassign the cell to Fibroblast (matching its UMAP coordinates)
seu_clean$cell_type[stray_epi_idx] <- "Fibroblast"
seu_clean$manual_annotation <- as.character(seu_clean$manual_annotation)
seu_clean$manual_annotation[stray_epi_idx] <- "Fibroblast"

# Re-apply factor levels to preserve plotting order and colors
seu_clean$cell_type <- factor(seu_clean$cell_type, levels = plot_order)
Idents(seu_clean) <- seu_clean$cell_type

# 4. Generate Abbreviations Based on the Corrected Cell Types
# Map abbreviations and use unname() to prevent Seurat metadata barcode mismatch errors
seu_clean$cell_type_abbr <- unname(cell_abbr[as.character(seu_clean$cell_type)])

# Re-map the color palette keys to use the new abbreviations
cell_cols_abbr <- setNames(cell_cols, cell_abbr[names(cell_cols)])

# 5. Define Mini-Axes and Themes for UMAPs
origin_x <- -15
origin_y <- -15
len      <- 5   

mini_axes <- data.frame(
  x    = c(origin_x, origin_x),
  y    = c(origin_y, origin_y),
  xend = c(origin_x + len, origin_x),
  yend = c(origin_y, origin_y + len)
)

# Common theme to apply to both UMAPs to keep code DRY
theme_umap <- theme(
  axis.line   = element_blank(),
  axis.ticks  = element_blank(),
  axis.text   = element_blank(),
  axis.title  = element_blank(),
  plot.margin = margin(t = 5, r = 5, b = -25, l = 30)
)

# 6. Create Combined UMAP (Full Names)
p_umap_combined <- DimPlot(
  seu_clean,
  reduction  = "umap.rpca",
  label      = TRUE,
  label.size = 5,
  cols       = cell_cols          
) +
  NoGrid() + NoLegend() + coord_cartesian(clip = "off") + theme_umap +
  geom_segment(
    data = mini_axes, aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE, arrow = arrow(type = "closed", length = unit(3, "pt")), linewidth = 0.4
  ) +
  annotate("text", x = origin_x + len, y = origin_y, label = "UMAP1", vjust = 1.5, size = 3) +
  annotate("text", x = origin_x - 0.35, y = origin_y + len, label = "UMAP2", hjust = 0.5, angle = 90, size = 3)

# 7. Create Split UMAP (Abbreviations & Clean Labels)
# Extract UMAP coordinates and merge with metadata to calculate manual label positions
umap_coords <- as.data.frame(Embeddings(seu_clean, reduction = "umap.rpca"))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
meta_df <- cbind(seu_clean@meta.data, umap_coords)

# Calculate label positions per condition, dropping groups with 0 cells
label_df <- meta_df %>%
  group_by(condition, cell_type_abbr) %>%
  summarize(
    n_cells = n(),
    UMAP1   = median(UMAP1),
    UMAP2   = median(UMAP2),
    .groups = "drop"
  ) %>%
  filter(n_cells > 0) # Safely drops groups with 0 cells

p_umap_split <- DimPlot(
  seu_clean,
  reduction  = "umap.rpca",
  split.by   = "condition",
  group.by   = "cell_type_abbr", 
  label      = FALSE,            # Turn off default Seurat labeling to prevent ghost labels
  cols       = cell_cols_abbr    
) +
  ggtitle(NULL) +                # Remove the "cell_type_abbr" global title
  geom_text(                     # Apply manual labels perfectly centered
    data = label_df, 
    aes(x = UMAP1, y = UMAP2, label = cell_type_abbr), 
    size = 5
  ) +
  NoGrid() + 
  NoLegend() + 
  coord_cartesian(clip = "off") + 
  theme_umap +
  geom_segment(
    data = mini_axes, aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE, arrow = arrow(type = "closed", length = unit(3, "pt")), linewidth = 0.4
  ) +
  annotate("text", x = origin_x + len, y = origin_y, label = "UMAP1", vjust = 1.5, size = 3) +
  annotate("text", x = origin_x - 0.35, y = origin_y + len, label = "UMAP2", hjust = 0.5, angle = 90, size = 3)

# 8. Calculate Relative Abundances for Bar Plot
# tab_counts <- table(CellType = Idents(seu_clean), Condition = seu_clean$condition)
# tab_prop <- prop.table(tab_counts, margin = 2)

# df_prop <- as.data.frame(tab_prop) %>% dplyr::rename(Proportion = Freq)

# Process proportions, calculate midpoints for labels, and format percentages
# df_bar_data <- df_prop %>%
  # mutate(
    # Percentage = round(Proportion * 100, 1),
    # Label = paste0(cell_abbr[as.character(CellType)], "=", Percentage, "%") 
  # ) %>%
  # group_by(Condition) %>%
  # arrange(Condition, desc(CellType)) %>% 
  # mutate(
    # cumulative = cumsum(Proportion),
    # midpoint = cumulative - (Proportion / 2)
  # ) %>%
  # ungroup()

# Generate staggered label positions to prevent text overlap in the bar chart
# df_labels <- df_bar_data %>% 
  # group_by(Condition) %>%
  # arrange(Condition, midpoint) %>% 
  # mutate(
    # label_pos = seq(0.0, 0.88, length.out = n()),
    # v_pos = case_when(
      # row_number() %% 3 == 1 ~ 0.65, 
      # row_number() %% 3 == 2 ~ 0.40, 
      # row_number() %% 3 == 0 ~ 0.15  
    # )
  # ) %>%
  # ungroup()
# 8. Calculate Sample-Level Proportions and Condition-Level Statistics
sample_meta <- seu_clean@meta.data %>%
  dplyr::select(sample, condition) %>%
  dplyr::distinct()

df_sample_props <- seu_clean@meta.data %>%
  dplyr::count(sample, cell_type) %>%
  tidyr::complete(sample, cell_type, fill = list(n = 0)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    total_cells = sum(n),
    proportion = n / total_cells,
    percentage = proportion * 100
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(sample_meta, by = "sample")

df_condition_stats <- df_sample_props %>%
  dplyr::group_by(condition, cell_type) %>%
  dplyr::summarize(
    mean_prop = mean(proportion),
    sd_prop   = sd(proportion),
    mean_pct  = mean(percentage),
    sd_pct    = sd(percentage),
    .groups   = "drop"
  )

# Save statistics for downstream formal testing
write.csv(df_condition_stats, "04-results/relative_abundance_stats.csv", row.names = FALSE)

# 9. Create Bar Plot
p_bar <- ggplot(df_bar_data, aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "stack", color = "black", width = 0.4) + 
  geom_point(data = df_labels, aes(x = v_pos, y = label_pos, color = CellType), size = 3) +
  geom_text(data = df_labels, aes(x = v_pos, y = label_pos, label = Label), hjust = 0, nudge_y = 0.015, size = 5) +
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

# 10. Combine Split Plot via Patchwork and Save Both Outputs
p_final <- p_umap_split / p_bar + plot_layout(heights = c(5, 1))

# Save the split panel exactly as before
ggsave(
  filename = "04-results/umap_rel_abundance.png",
  plot     = p_final,
  width    = 24,    
  height   = 10,
  units    = "in",
  dpi      = 300
)

# Save the combined UMAP independently, mathematically matched to a single split panel
ggsave(
  filename = "04-results/umap_combined.png",
  plot     = p_umap_combined,
  width    = 8,       # 24 inches total width / 3 panels = 8 inches wide
  height   = 8.33,    # 10 inches total height * (5/6 ratio for top panel) = 8.33 inches high
  units    = "in",
  dpi      = 300
)

saveRDS(seu_clean, "03-analysis_scratch/seu_after_visualization.rds")


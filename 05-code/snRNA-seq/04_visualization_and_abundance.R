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
library(ggtext)

# 1. Load Data & Ensure Ordering ----
seu_clean <- readRDS("03-analysis_scratch/seu_annotated.rds")

# Lock the condition levels so the facets and bars are ordered correctly
seu_clean$condition <- factor(
  seu_clean$condition,
  levels = c("Control", "MCT-Water", "MCT-Blumeria")
)

# Initialize standard cell types from manual annotations
Idents(seu_clean) <- "manual_annotation"
seu_clean$cell_type <- Idents(seu_clean)

# 2. Define Plot Aesthetics (Order, Colors, Abbreviations) ----
celltype_order <- c(
  "Cardiomyocyte", "Fibroblast", "Macrophage", "Monocyte",
  "Dendritic Cell", "Neutrophil", "T Cell", "B Cell", "NK Cell",
  "Pericyte", "Vascular EC", "Endocardial EC", "Lymphatic EC",
  "Neuronal"
)

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

# Reverse standard order so items stack predictably in ggplot bars
plot_order <- rev(celltype_order)  

# 3. Apply Clean Factor Levels ----
seu_clean$cell_type <- factor(seu_clean$cell_type, levels = plot_order)
Idents(seu_clean) <- seu_clean$cell_type

# 4. Generate Abbreviations Based on Cell Types
# Map abbreviations and use unname() to prevent Seurat metadata barcode mismatch errors
seu_clean$cell_type_abbr <- unname(cell_abbr[as.character(seu_clean$cell_type)])

# Re-map the color palette keys to use the new abbreviations
cell_cols_abbr <- setNames(cell_cols, cell_abbr[names(cell_cols)])

# 5. Define Mini-Axes and Themes for UMAPs ----
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

# 6. Create Combined UMAP (Full Names) ----
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

# 7. Create Split UMAP (Abbreviations & Clean Labels) ----
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

# Define clean display labels using markdown tags
condition_labels <- c(
  "Control"      = "Control",
  "MCT-Water"    = "MCT-Water",
  "MCT-Blumeria" = "MCT-*Blumeria*" # Or "MCT-*Blumeria*"
)

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
  facet_wrap(
    ~ condition, 
    nrow = 1, 
    labeller = as_labeller(condition_labels)
  ) +
  theme(strip.text = element_markdown(face = "bold", size = 14))+
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

# 8. Calculate Sample-Level Proportions and Condition-Level Statistics ----
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

# 9. Format Data for Bar Plot ----
df_bar_data <- df_condition_stats %>%
  # Rename columns to seamlessly fit your existing ggplot code
  dplyr::rename(
    Condition  = condition, 
    CellType   = cell_type, 
    Proportion = mean_prop, 
    Percentage = mean_pct
  ) %>%
  group_by(Condition) %>%
  # Normalize mean proportions so the stacked bars reach exactly 100% 
  # (averaging replicates can sometimes equal 0.999 or 1.001)
  mutate(Proportion = Proportion / sum(Proportion)) %>% 
  arrange(Condition, desc(CellType)) %>% 
  mutate(
    Percentage = round(Proportion * 100, 1),
    Label = paste0(cell_abbr[as.character(CellType)], "=", Percentage, "%")#,
    # cumulative = cumsum(Proportion),
    # midpoint = cumulative - (Proportion / 2)
  ) %>%
  ungroup()

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
# # Map cell types into a strict 3-Column x 5-Row Matrix below each bar
# df_labels <- df_bar_data %>% 
#   group_by(Condition) %>%
#   # Ensure items follow celltype_order top-to-bottom, left-to-right
#   mutate(type_order = match(as.character(CellType), celltype_order)) %>%
#   arrange(Condition, type_order) %>%
#   mutate(
#     idx     = row_number(),
#     # 3 Columns: Col 1 (1..5), Col 2 (6..10), Col 3 (11..14)
#     col_num = ((idx - 1) %/% 5) + 1,
#     row_num = ((idx - 1) %% 5) + 1,
#     
#     # Fixed Column Coordinates (Horizontal spacing across facet)
#     grid_col_point = case_when(
#       col_num == 1 ~ 0.01,
#       col_num == 2 ~ 0.35,
#       col_num == 3 ~ 0.68
#     ),
#     grid_col_text = grid_col_point + 0.035,
#     
#     # Fixed Row Coordinates (Vertical spacing below the bar)
#     grid_row = case_when(
#       row_num == 1 ~ 0.65,
#       row_num == 2 ~ 0.48,
#       row_num == 3 ~ 0.31,
#       row_num == 4 ~ 0.14,
#       row_num == 5 ~ -0.03
#     )
#   ) %>%
#   ungroup()

# Map 14 cell types into a 4-Column Grid (Filled ROW BY ROW)
# Row 1: CM, FB, MФ, Mono
# Row 2: DC, Neu, T, B
# Row 3: NK, Peri, vEC, eEC
# Row 4: lEC, Neur
df_labels <- df_bar_data %>% 
  group_by(Condition) %>%
  mutate(type_order = match(as.character(CellType), celltype_order)) %>%
  arrange(Condition, type_order) %>%
  mutate(
    idx     = row_number(),
    col_num = ((idx - 1) %% 4) + 1,
    row_num = ((idx - 1) %/% 4) + 1,
    
    # 4 Column Starting Points (evenly spaced at 0%, 25%, 50%, and 75% width)
    grid_col_point = case_when(
      col_num == 1 ~ 0.01,
      col_num == 2 ~ 0.26,
      col_num == 3 ~ 0.51,
      col_num == 4 ~ 0.76
    ),
    grid_col_text = grid_col_point + 0.032,
    
    # 4 Fixed Row Heights below the bar
    grid_row = case_when(
      row_num == 1 ~ 0.60,
      row_num == 2 ~ 0.42,
      row_num == 3 ~ 0.24,
      row_num == 4 ~ 0.06
    )
  ) %>%
  ungroup()

# 10. Create Bar Plot ----
# p_bar <- ggplot(df_bar_data, aes(x = 1, y = Proportion, fill = CellType)) +
  # geom_col(position = "stack", color = "black", width = 0.4) + 
  # geom_point(data = df_labels, aes(x = v_pos, y = label_pos, color = CellType), size = 3) +
  # geom_text(data = df_labels, aes(x = v_pos, y = label_pos, label = Label), hjust = 0, nudge_y = 0.015, size = 5) +
  # coord_flip(xlim = c(0.2, 1.2), ylim = c(0, 1), clip = "off") + 
  # facet_wrap(~ Condition, nrow = 1) +
  # scale_fill_manual(values = cell_cols, guide = "none") +
  # scale_color_manual(values = cell_cols, guide = "none") +
  # labs(title = "Relative Abundance (%)", x = NULL, y = NULL) +
  # theme_minimal() +
  # theme(
  #   plot.title       = element_text(hjust = 0.5),
  #   strip.background = element_blank(),
  #   strip.text       = element_blank(),
  #   panel.background = element_blank(),
  #   plot.background  = element_blank(),
  #   panel.grid       = element_blank(),
  #   axis.text.x      = element_blank(),
  #   axis.ticks.x     = element_blank(),
  #   axis.text.y      = element_blank(),
  #   axis.ticks.y     = element_blank()
  # )

p_bar <- ggplot(df_bar_data, aes(x = 1, y = Proportion, fill = CellType)) +
  geom_col(position = "stack", color = "black", width = 0.35) + 
  
  # Colored dots at exact grid coordinates
  geom_point(
    data = df_labels, 
    aes(x = grid_row, y = grid_col_point, color = CellType), 
    size = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Text labels at exact grid coordinates
  geom_text(
    data = df_labels, 
    aes(x = grid_row, y = grid_col_text, label = Label), 
    hjust = 0, 
    vjust = 0.5,
    size = 3.8,
    inherit.aes = FALSE
  ) +
  
  # Coordinate flip with expanded xlim to fit the 5-row grid cleanly
  coord_flip(xlim = c(-0.12, 1.25), ylim = c(0, 1), clip = "off") + 
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  scale_color_manual(values = cell_cols, guide = "none") +
  labs(title = "Relative Abundance (%)", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    plot.margin      = margin(t = 5, r = 10, b = 10, l = 10)
  )

# 11. Combine Figures and Export Outputs ----
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

# 12. Supplementary Nuclei Summary Table ----
sample_order <- c(
  "Control11", "Control12", "Control14", "Control15",
  "Water1", "Water5", "Water6", "Water9",
  "Blumeria1", "Blumeria3", "Blumeria4", "Blumeria7"
)

# 1. Summarize and force column order
celltype_sample_summary <- seu_clean@meta.data %>%
  group_by(cell_type, sample) %>%
  summarize(n_nuclei = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = n_nuclei, values_fill = 0) %>%
  select(cell_type, all_of(sample_order)) # Ensures correct left-to-right sample ordering

celltype_totals <- seu_clean@meta.data %>%
  group_by(cell_type) %>%
  summarize(
    Total_Nuclei = n(),
    Percent_of_Total = round((n() / nrow(seu_clean@meta.data)) * 100, 2),
    .groups = "drop"
  )

master_table <- left_join(celltype_totals, celltype_sample_summary, by = "cell_type") %>%
  mutate(cell_type = as.character(cell_type))

# 2. Build total summary row
sample_cols <- setdiff(colnames(master_table), c("cell_type", "Total_Nuclei", "Percent_of_Total"))

total_row <- data.frame(
  cell_type = "Total",
  Total_Nuclei = sum(master_table$Total_Nuclei),
  Percent_of_Total = round(sum(master_table$Percent_of_Total), 2)
)

for (col in sample_cols) {
  total_row[[col]] <- sum(master_table[[col]])
}

# 3. Combine and export cleanly ordered table
master_table_final <- bind_rows(master_table, total_row)
write.csv(master_table_final, "04-results/Supplementary_Table_Nuclei_Counts_per_CellType.csv", row.names = FALSE)

# Export updated seurat object ----
saveRDS(seu_clean, "03-analysis_scratch/seu_after_visualization.rds")

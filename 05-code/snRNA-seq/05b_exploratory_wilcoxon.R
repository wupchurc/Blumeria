# ==============================================================================
# Script: 05b_exploratory_wilcoxon.R
# Purpose: Exploratory single-cell DEG analysis (Wilcoxon) and MA jitter plots.
#          NOTE: For final manuscript statistics, refer to 05_pseudobulk_DESeq2.R
# Input:   03-analysis_scratch/seu_after_visualization.rds
# Output:  03-analysis_scratch/DEG_Wilcox_*.rds
# ==============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)

# 1. Load updated Seurat object
seu_obj <- readRDS("03-analysis_scratch/seu_after_visualization.rds")

# 3. Single-Cell Wilcoxon DEG Function
run_wilcox_deg <- function(seu_obj, cell_type, min_pct = 0.1, logfc_thresh = 0, save_results = FALSE) {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) stop(paste("No cells found for cell type:", cell_type))
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n")
  
  Idents(seu_subset) <- "condition"
  
  # Helper to run FindMarkers and standardize to DESeq2 naming conventions
  get_contrast <- function(ident_1, ident_2) {
    res <- FindMarkers(
      seu_subset, 
      ident.1 = ident_1, 
      ident.2 = ident_2, 
      test.use = "wilcox",
      logfc.threshold = logfc_thresh, 
      min.pct = min_pct               
    )
    
    res <- res %>%
      rownames_to_column(var = "gene") %>%
      dplyr::rename(log2FoldChange = avg_log2FC, padj = p_val_adj, pvalue = p_val) %>%
      column_to_rownames(var = "gene")
    
    return(res)
  }
  
  cat("Running Wilcoxon contrasts...\n")
  res_water_vs_ctrl <- get_contrast("MCT-Water", "Control")
  res_blum_vs_ctrl  <- get_contrast("MCT-Blumeria", "Control")
  res_blum_vs_water <- get_contrast("MCT-Blumeria", "MCT-Water")
  res_water_vs_blum <- get_contrast("MCT-Water", "MCT-Blumeria")
  
  results_list <- list(
    water_vs_ctrl = res_water_vs_ctrl,
    blum_vs_ctrl  = res_blum_vs_ctrl,
    blum_vs_water = res_blum_vs_water,
    water_vs_blum = res_water_vs_blum
  )
  
  if (save_results) {
    filename <- paste0("DEG_Wilcox_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = paste0("03-analysis_scratch/", filename))
  }
  
  return(results_list)
}

# 4. Execute Wilcoxon Across All Cell Types
deg_results <- list()
cell_types <- levels(seu_obj)  
for (cell_type in cell_types) {
  cat("\n=== Processing", cell_type, "===\n")
  deg_results[[cell_type]] <- run_wilcox_deg(seu_obj, cell_type, save_results = TRUE)
}

# 5. Plotting Function (Refined)
create_ma_plot <- function(deg_results, cell_types, contrast_name, padj_thresh = 0.05, lfc_thresh = 0.5) {
  
  plot_data <- data.frame()
  up_counts <- integer(length(cell_types))
  down_counts <- integer(length(cell_types))
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    if (contrast_name == "MCT-Water vs Control") {
      res <- deg_results[[cell_type]]$water_vs_ctrl
    } else if (contrast_name == "MCT-Water vs MCT-Blumeria") {
      res <- deg_results[[cell_type]]$water_vs_blum
    } else if (contrast_name == "MCT-Blumeria vs Control") {
      res <- deg_results[[cell_type]]$blum_vs_ctrl
    }
    
    df <- as.data.frame(res) %>%
      mutate(
        gene = rownames(.),
        cell_type = cell_type,
        cell_type_num = i,  
        significant = !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) > lfc_thresh,
        direction = case_when(
          significant & log2FoldChange > 0 ~ "UP",
          significant & log2FoldChange < 0 ~ "DOWN",
          TRUE ~ "NS"
        )
      )
    
    plot_data <- rbind(plot_data, df)
    up_counts[i] <- sum(df$direction == "UP", na.rm = TRUE)
    down_counts[i] <- sum(df$direction == "DOWN", na.rm = TRUE)
  }
  
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types)
  display_title <- gsub("Blumeria", "*Blumeria*", contrast_name)
  
  y_max <- 5
  y_min <- -5
  
  p <- ggplot(plot_data, aes(x = cell_type_num, y = log2FoldChange)) +
    coord_cartesian(ylim = c(y_min, y_max)) + # Keeps visual range without deleting outliers
    geom_jitter(data = filter(plot_data, direction == "NS"), width = 0.25, height = 0, size = 0.3, alpha = 0.2, color = "gray70") +
    geom_jitter(data = filter(plot_data, direction == "UP"), width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#E31A1C") +
    geom_jitter(data = filter(plot_data, direction == "DOWN"), width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#1F78B4") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_continuous(breaks = 1:length(cell_types), labels = cell_types) +
    labs(title = display_title, x = NULL, y = "Log2FoldChange") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 16) 
    )
  
  for (i in seq_along(cell_types)) {
    p <- p + annotate("text", x = i, y = y_max*0.95, label = up_counts[i], color = "#E31A1C", size = 3, fontface = "bold")
    p <- p + annotate("text", x = i, y = y_min*0.95, label = down_counts[i], color = "#1F78B4", size = 3, fontface = "bold")
  }
  
  p <- p + annotate("text", x = length(cell_types) + 0.8, y = y_max * 0.95, label = "UP", color = "#E31A1C", size = 4, fontface = "bold") +
    annotate("text", x = length(cell_types) + 0.8, y = y_min * 0.95, label = "DOWN", color = "#1F78B4", size = 4, fontface = "bold")
  
  return(p)
}

# 6. Generate Plots
p_water_blum <- create_ma_plot(deg_results, rev(cell_types), "MCT-Water vs MCT-Blumeria", padj_thresh = 0.05)
p_water_ctrl <- create_ma_plot(deg_results, rev(cell_types), "MCT-Water vs Control", padj_thresh = 0.05)
p_combined <- (p_water_ctrl | p_water_blum)

# Setup layout margins and turn off clipping so right-side text (UP/DOWN) renders outside the bounds
p_combined <- p_combined +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, r = 30, b = 25, l = 10))

# Export combined and individual plots
ggsave("04-results/Differential_Expression_Plot_Wilcox.png", plot = p_combined, width = 16, height = 7, units = "in", dpi = 300)

# ==============================================================================
# 7. Export Targeted ShinyGO Gene Lists (MCT-Water vs MCT-Blumeria)
# ==============================================================================

# Helper function to extract gene symbols and write .txt files
export_shinygo_lists <- function(res, cell_name, comp_name, lfc_thresh = 0.5, padj_thresh = 0.05) {
  if (is.null(res)) return(NULL)
  
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  # Filter significant features
  sig_up   <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange > lfc_thresh]
                  sig_down <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange < -lfc_thresh]
                                  background <- df$gene
                                  
                                  # Create directory for output
                                  out_dir <- "03-analysis_scratch/shinygo_wilcox"
                                  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                                  
                                  # Write text files
                                  writeLines(sig_up,     sprintf("%s/%s_%s_up_genes.txt", out_dir, cell_name, comp_name))
                                  writeLines(sig_down,   sprintf("%s/%s_%s_down_genes.txt", out_dir, cell_name, comp_name))
                                  writeLines(background, sprintf("%s/%s_%s_background.txt", out_dir, cell_name, comp_name))
                                  
                                  cat(sprintf("Exported %s (%s): %d UP, %d DOWN\n", 
                                              cell_name, comp_name, length(sig_up), length(sig_down)))
}

# Define target cell types
target_cells <- c("Cardiomyocyte", "Fibroblast", "Macrophage")

cat("\n=== Exporting Targeted ShinyGO Lists (water_vs_blum) ===\n")

for (cell_type in target_cells) {
  if (cell_type %in% names(deg_results)) {
    clean_cell_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
    res_water_vs_blum <- deg_results[[cell_type]]$water_vs_blum
    
    export_shinygo_lists(res_water_vs_blum, clean_cell_name, "water_vs_blum")
  } else {
    warning(sprintf("Cell type '%s' not found in deg_results!", cell_type))
  }
}

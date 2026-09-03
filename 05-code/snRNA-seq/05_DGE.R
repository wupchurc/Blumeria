# ==============================================================================
# Script: 05_  .R
# Purpose: 
# Input:   03-analysis_scratch/seu_after_visualization.rds
# Output:  
# ==============================================================================

# load libraries 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)
library(EnhancedVolcano)
library(ggtext)

# Load the updated Seurat object from the previous step
seu_obj <- readRDS("03-analysis_scratch/seu_after_visualization.rds")

# Ensure the active identity is set to your cleaned cell types
Idents(seu_obj) <- seu_obj$cell_type

# ---- Consolidated Pseudobulking and DESeq2 Function ----
run_pseudobulk_deg <- function(seu_obj, cell_type, min_counts = 10, alpha = 0.05,
                               plot_title = NULL, save_results = FALSE) {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) {
    stop(paste("No cells found for cell type:", cell_type))
  }
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n")
  
  cts <- AggregateExpression(seu_subset, assays = "RNA", group.by = "orig.ident",
                             slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  colData <- seu_subset@meta.data %>%
    dplyr::select(orig.ident, condition) %>%
    dplyr::distinct() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "orig.ident")
  
  colData <- colData[colnames(cts), , drop = FALSE]
  
  # Ensure condition is a factor and drop any empty levels (like Blumeria for Epi cells)
  colData$condition <- droplevels(factor(colData$condition))
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = colData,
    design = ~ condition)
  
  keep <- rowSums(counts(dds)) >= min_counts
  dds <- dds[keep, ]
  cat("Filtered to", nrow(dds), "genes (min total counts =", min_counts, ")\n")
  
  # Only relevel if Control is actually present in this cell type
  if ("Control" %in% levels(dds$condition)) {
    dds$condition <- relevel(dds$condition, ref = "Control")
  }
  
  dds <- DESeq(dds)
  
  # Helper function to safely extract results only if both conditions exist
  safe_results <- function(c1, c2) {
    if (c1 %in% levels(dds$condition) & c2 %in% levels(dds$condition)) {
      return(results(dds, contrast = c("condition", c1, c2), alpha = alpha))
    } else {
      cat(sprintf("  -> Skipping %s vs %s (missing cells in one condition)\n", c1, c2))
      return(NULL)
    }
  }
  
  res_water_vs_ctrl <- safe_results("MCT-Water", "Control")
  res_blum_vs_ctrl  <- safe_results("MCT-Blumeria", "Control")
  res_blum_vs_water <- safe_results("MCT-Blumeria", "MCT-Water")
  res_water_vs_blum <- safe_results("MCT-Water", "MCT-Blumeria")
  
  # Save results if requested
  if (save_results) {
    results_list <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl  = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum
    )
    filename <- paste0("DEG_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = paste0("03-analysis_scratch/", filename))
  }
  
  return(
    list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl  = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum,
      dds = dds
    )
  )
}

# ---- MA Plots of up and down regulated gene counts ----

# Run DESeq2 on all cell types
deg_results <- list()
cell_types <- levels(seu_obj)  
for (cell_type in cell_types) {
  cat("\n=== Processing", cell_type, "===\n")
  deg_results[[cell_type]] <- run_pseudobulk_deg(seu_obj, cell_type, 
                                                 alpha = 0.2, save_results = FALSE)
}

# Function to create one panel of MA plots
create_ma_plot <- function(deg_results, cell_types, contrast_name, 
                           padj_thresh = 0.05, lfc_thresh = 0.5) {
  
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
    
    # NEW: Skip this iteration if the contrast wasn't possible for this cell type
    if (is.null(res)) { 
      up_counts[i] <- 0
      down_counts[i] <- 0
      next 
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
  display_title <- gsub("Blumeria", "<i>Blumeria</i>", contrast_name)
  
  p <- ggplot(plot_data, aes(x = cell_type_num, y = log2FoldChange)) +
    ylim(-10,10)+
    geom_jitter(data = filter(plot_data, direction == "NS"),
                width = 0.25, height = 0, size = 0.3, alpha = 0.2, color = "gray70") +
    geom_jitter(data = filter(plot_data, direction == "UP"),
                width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#E31A1C") +
    geom_jitter(data = filter(plot_data, direction == "DOWN"),
                width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#1F78B4") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_continuous(breaks = 1:length(cell_types),
                       labels = cell_types) +
    labs(title = display_title, 
         x = NULL,
         y = "Log2FoldChange") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 16) 
    )
  
  y_max <- 10
  y_min <- -10
  
  for (i in seq_along(cell_types)) {
    p <- p + annotate("text", x = i, y = y_max*0.95,
                      label = up_counts[i], color = "#E31A1C", size = 3, fontface = "bold")
    p <- p + annotate("text", x = i, y = y_min*0.95,
                      label = down_counts[i], color = "#1F78B4", size = 3, fontface = "bold")
  }
  
  p <- p + annotate("text", x = length(cell_types) + 0.8, y = y_max * 0.95,
                    label = "UP", color = "#E31A1C", size = 4, fontface = "bold") +
    annotate("text", x = length(cell_types) + 0.8, y = y_min * 0.95,
             label = "DOWN", color = "#1F78B4", size = 4, fontface = "bold")
  
  return(p)
}

# Create plots for each contrast
p_water_ctrl <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Water vs Control")
p_blum_ctrl <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Blumeria vs Control")
p_water_blum <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Water vs MCT-Blumeria")

# Combine using patchwork
p_combined_paper <- (p_water_ctrl | p_blum_ctrl)

print(p_combined_paper)

ggsave(
  filename = "04-results/Differential_Expression_Plot_p_05.png",
  plot = p_combined_paper,
  width = 12,
  height = 8,
  dpi = 300
)

p <- p_combined_paper +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, r = 30, b = 25, l = 10))

ggsave("04-results/Differential_Expression_Plot.png", plot = p, width = 16, height = 7, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_water_vs_ctrl.png", plot = p_water_ctrl, width = 7, height = 4, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_blum_vs_ctrl.png", plot = p_blum_ctrl, width = 7, height = 4, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_water_vs_blum.png", plot = p_water_blum, width = 7, height = 4, units = "in", dpi = 300)

# ----
run_pseudobulk_deg(seu_obj, "Cardiomyocyte", alpha = 0.2, save_results = TRUE)
run_pseudobulk_deg(seu_obj, "Macrophage", alpha = 0.2, save_results = TRUE)

# ---- Export Significant Genes for ShinyGO ----
export_shinygo_lists <- function(res, cell_name, comp_name, lfc_thresh = 0.5, padj_thresh = 0.05) {
  if (is.null(res)) return(NULL)
  
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  sig_up <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange > lfc_thresh]
  sig_down <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange < -lfc_thresh]
  background <- df$gene
  
  writeLines(sig_up, sprintf("03-analysis_scratch/%s_%s_up_genes.txt", cell_name, comp_name))
  writeLines(sig_down, sprintf("03-analysis_scratch/%s_%s_down_genes.txt", cell_name, comp_name))
  writeLines(background, sprintf("03-analysis_scratch/%s_%s_background.txt", cell_name, comp_name))
}

# Reload the saved RDS files if they aren't in your environment
cm_results <- readRDS("03-analysis_scratch/DEG_Cardiomyocyte.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophage.rds")

# Export lists for Cardiomyocytes
export_shinygo_lists(cm_results$water_vs_ctrl, "cm", "water_vs_ctrl")
export_shinygo_lists(cm_results$water_vs_blum, "cm", "water_vs_blum")
export_shinygo_lists(cm_results$blum_vs_ctrl, "cm", "blum_vs_ctrl")

# Export lists for Macrophages
export_shinygo_lists(mac_results$water_vs_ctrl, "mac", "water_vs_ctrl")
export_shinygo_lists(mac_results$water_vs_blum, "mac", "water_vs_blum")
export_shinygo_lists(mac_results$blum_vs_ctrl, "mac", "blum_vs_ctrl")
# ==============================================================================
# Script: 05_DGE.R
# Purpose: Perform pseudobulk differential gene expression (DGE) analysis across 
#          all cell types using DESeq2. Generate comparative MA plots and export 
#          significant gene lists for downstream pathway enrichment (e.g., ShinyGO).
# Input:   03-analysis_scratch/seu_after_visualization.rds
# Output:  04-results/Differential_Expression_Plot*.png
#          03-analysis_scratch/DEG_*.rds (DESeq2 results objects)
#          03-analysis_scratch/*_genes.txt (Gene lists for ShinyGO)
# ==============================================================================

# ---- Load Required Libraries ---- 
library(Seurat)
library(DESeq2)
library(tidyverse)
library(patchwork)
library(ggtext)

# ---- Load and Prepare Data ----
# Load the updated Seurat object from the visualization step
seu_obj <- readRDS("03-analysis_scratch/seu_after_visualization.rds")

# Ensure the active identity is set to your cleaned cell types
Idents(seu_obj) <- seu_obj$cell_type

# ---- Define Pseudobulking and DGE Function ----
# This function aggregates counts per sample for a given cell type and runs DESeq2
run_pseudobulk_deg <- function(seu_obj, cell_type, min_counts = 10, 
                               alpha = 0.05, save_results = FALSE) {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  # 1. Subset to the specific cell type
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) {
    stop(paste("No cells found for cell type:", cell_type))
  }
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n")
  
  # 2. Aggregate counts across cells to create a pseudobulk matrix by sample
  cts <- AggregateExpression(seu_subset, assays = "RNA", group.by = "orig.ident",
                             slot = "counts", return.seurat = FALSE)
  cts <- cts$RNA
  
  # 3. Extract and format sample metadata directly from the Seurat object
  colData <- seu_subset@meta.data %>%
    dplyr::select(orig.ident, condition) %>%
    dplyr::distinct() %>%
    tibble::remove_rownames() %>% # Remove cell barcodes
    tibble::column_to_rownames(var = "orig.ident")
  
  # Ensure the metadata rows perfectly align with the count matrix columns
  colData <- colData[colnames(cts), , drop = FALSE]
  
  # 4. Clean factor levels to prevent DESeq2 errors (e.g., if a condition has 0 cells)
  colData$condition <- droplevels(factor(colData$condition))
  
  # 5. Initialize DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = colData,
    design = ~ condition)
  
  # Filter out low-count genes to improve statistical power
  keep <- rowSums(counts(dds)) >= min_counts
  dds <- dds[keep, ]
  cat("Filtered to", nrow(dds), "genes (min total counts =", min_counts, ")\n")
  
  # 6. Set reference level to "Control" if it exists in this cell type
  if ("Control" %in% levels(dds$condition)) {
    dds$condition <- relevel(dds$condition, ref = "Control")
  }
  
  # 7. Run standard DESeq2 pipeline
  dds <- DESeq(dds)
  
  # 8. Helper function to extract contrasts safely
  # Avoids crashing if a specific condition group was dropped due to missing cells
  safe_results <- function(c1, c2) {
    if (c1 %in% levels(dds$condition) & c2 %in% levels(dds$condition)) {
      return(results(dds, contrast = c("condition", c1, c2), alpha = alpha))
    } else {
      cat(sprintf("  -> Skipping %s vs %s (missing cells in one condition)\n", c1, c2))
      return(NULL)
    }
  }
  
  # Extract all relevant pair-wise contrasts
  res_water_vs_ctrl <- safe_results("MCT-Water", "Control")
  res_blum_vs_ctrl  <- safe_results("MCT-Blumeria", "Control")
  res_blum_vs_water <- safe_results("MCT-Blumeria", "MCT-Water")
  res_water_vs_blum <- safe_results("MCT-Water", "MCT-Blumeria")
  
  # 9. Save raw results to disk if requested
  if (save_results) {
    results_list <- list(
      water_vs_ctrl = res_water_vs_ctrl,
      blum_vs_ctrl  = res_blum_vs_ctrl,
      blum_vs_water = res_blum_vs_water,
      water_vs_blum = res_water_vs_blum
    )
    # Sanitize cell type name for the filename (e.g., spaces to underscores)
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

# ---- Global Differential Expression Iteration ----
# Loop through all cell types and store temporary DESeq2 results for plotting
deg_results <- list()
cell_types <- levels(seu_obj)  
for (cell_type in cell_types) {
  cat("\n=== Processing", cell_type, "===\n")
  deg_results[[cell_type]] <- run_pseudobulk_deg(seu_obj, cell_type, 
                                                 alpha = 0.2, save_results = FALSE)
}

# ---- Plotting Function: Summary MA Plot ----
# Generates a customized 1D scatter plot (strip chart) of Log2FoldChanges
create_ma_plot <- function(deg_results, cell_types, contrast_name, 
                           padj_thresh = 0.05, lfc_thresh = 0.5) {
  
  plot_data <- data.frame()
  up_counts <- integer(length(cell_types))
  down_counts <- integer(length(cell_types))
  
  # Compile data across all cell types for the specified contrast
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    # Retrieve the correct results object based on the contrast requested
    if (contrast_name == "MCT-Water vs Control") {
      res <- deg_results[[cell_type]]$water_vs_ctrl
    } else if (contrast_name == "MCT-Water vs MCT-Blumeria") {
      res <- deg_results[[cell_type]]$water_vs_blum
    } else if (contrast_name == "MCT-Blumeria vs Control") {
      res <- deg_results[[cell_type]]$blum_vs_ctrl
    }
    
    # Skip iteration entirely if this contrast couldn't be calculated
    if (is.null(res)) { 
      up_counts[i] <- 0
      down_counts[i] <- 0
      next 
    }
    
    # Format results and assign significance categories
    df <- as.data.frame(res) %>%
      mutate(
        gene = rownames(.),
        cell_type = cell_type,
        cell_type_num = i,  
        significant = !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) > lfc_thresh,
        direction = case_when(
          significant & log2FoldChange > 0 ~ "UP",
          significant & log2FoldChange < 0 ~ "DOWN",
          TRUE ~ "NS" # Not significant
        )
      )
    
    plot_data <- rbind(plot_data, df)
    
    # Tally up/down genes for annotation
    up_counts[i] <- sum(df$direction == "UP", na.rm = TRUE)
    down_counts[i] <- sum(df$direction == "DOWN", na.rm = TRUE)
  }
  
  # Ensure cell types remain in the specified order on the X-axis
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types)
  # Italicize Blumeria using ggtext markdown
  display_title <- gsub("Blumeria", "<i>Blumeria</i>", contrast_name)
  
  # Build the ggplot
  p <- ggplot(plot_data, aes(x = cell_type_num, y = log2FoldChange)) +
    ylim(-10,10) +
    # Layer points: NS first, then significant points on top
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
  
  # Add numerical count annotations at the top and bottom of the plot
  y_max <- 10
  y_min <- -10
  
  for (i in seq_along(cell_types)) {
    p <- p + annotate("text", x = i, y = y_max*0.95,
                      label = up_counts[i], color = "#E31A1C", size = 3, fontface = "bold")
    p <- p + annotate("text", x = i, y = y_min*0.95,
                      label = down_counts[i], color = "#1F78B4", size = 3, fontface = "bold")
  }
  
  # Add row labels (UP/DOWN) to the far right margin (+1.2 offset to prevent crowding)
  p <- p + annotate("text", x = length(cell_types) + 0.8, y = y_max * 0.95,
                    label = "UP", color = "#E31A1C", size = 4, fontface = "bold") +
    annotate("text", x = length(cell_types) + 0.8, y = y_min * 0.95,
             label = "DOWN", color = "#1F78B4", size = 4, fontface = "bold")
  
  return(p)
}

# ---- Generate and Save Plots ----
# Create plots for each contrast
p_water_ctrl <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Water vs Control")
p_blum_ctrl <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Blumeria vs Control")
p_water_blum <- create_ma_plot(deg_results, rev(cell_types), padj_thresh = 0.05, lfc_thresh = 0.5, "MCT-Water vs MCT-Blumeria")

# Combine main comparators side-by-side
p_combined <- (p_water_ctrl | p_blum_ctrl)

# Setup layout margins and turn off clipping so right-side text (UP/DOWN) renders outside the bounds
p_combined <- p_combined +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, r = 30, b = 25, l = 10))

# Export combined and individual plots
ggsave("04-results/Differential_Expression_Plot.png", plot = p_combined, width = 16, height = 7, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_water_vs_ctrl.png", plot = p_water_ctrl, width = 7, height = 4, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_blum_vs_ctrl.png", plot = p_blum_ctrl, width = 7, height = 4, units = "in", dpi = 300)
ggsave("04-results/Differential_Expression_Plot_water_vs_blum.png", plot = p_water_blum, width = 7, height = 4, units = "in", dpi = 300)

# ---- Run and Save Target Cell Type Results ----
# Run full analysis on target clusters and save RDS files for downstream processing
run_pseudobulk_deg(seu_obj, "Cardiomyocyte", alpha = 0.2, save_results = TRUE)
run_pseudobulk_deg(seu_obj, "Macrophage", alpha = 0.2, save_results = TRUE)

# ---- Export Significant Genes for ShinyGO ----
# Helper function to extract gene symbols based on thresholds and write .txt files
export_shinygo_lists <- function(res, cell_name, comp_name, lfc_thresh = 0.5, padj_thresh = 0.05) {
  if (is.null(res)) return(NULL)
  
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  # Filter significant features
  sig_up <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange > lfc_thresh]
  sig_down <- df$gene[!is.na(df$padj) & df$padj < padj_thresh & df$log2FoldChange < -lfc_thresh]
  background <- df$gene
  
  # Write raw text files
  writeLines(sig_up, sprintf("03-analysis_scratch/%s_%s_up_genes.txt", cell_name, comp_name))
  writeLines(sig_down, sprintf("03-analysis_scratch/%s_%s_down_genes.txt", cell_name, comp_name))
  writeLines(background, sprintf("03-analysis_scratch/%s_%s_background.txt", cell_name, comp_name))
}

# Reload the saved RDS files generated by run_pseudobulk_deg
cm_results <- readRDS("03-analysis_scratch/DEG_Cardiomyocyte.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophage.rds")

# Execute exports for Cardiomyocytes
export_shinygo_lists(cm_results$water_vs_ctrl, "cm", "water_vs_ctrl")
export_shinygo_lists(cm_results$water_vs_blum, "cm", "water_vs_blum")
export_shinygo_lists(cm_results$blum_vs_ctrl, "cm", "blum_vs_ctrl")

# Execute exports for Macrophages
export_shinygo_lists(mac_results$water_vs_ctrl, "mac", "water_vs_ctrl")
export_shinygo_lists(mac_results$water_vs_blum, "mac", "water_vs_blum")
export_shinygo_lists(mac_results$blum_vs_ctrl, "mac", "blum_vs_ctrl")
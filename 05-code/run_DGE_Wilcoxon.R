# load libraries 
library(Seurat)
library(tidyverse)
library(patchwork)
library(EnhancedVolcano)
library(ggtext)

seu_obj <- readRDS("03-analysis_scratch/seu_for_DGE.rds")

# ---- NEW: Add Condition to Global Metadata ----
# Doing this once here makes the subsetting and grouping much easier
seu_obj$condition <- case_when(
  str_detect(seu_obj$orig.ident, "Blumeria") ~ "MCT-Blumeria",
  str_detect(seu_obj$orig.ident, "Water") ~ "MCT-Water",
  str_detect(seu_obj$orig.ident, "Control") ~ "Control"
)

# ---- Single-Cell Wilcoxon DEG Function ----
run_wilcox_deg <- function(seu_obj, cell_type, min_pct = 0.1, logfc_thresh = 0,
                           save_results = FALSE) {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  # Subset to specified cell type
  cells <- WhichCells(seu_obj, idents = cell_type)
  if (length(cells) == 0) {
    stop(paste("No cells found for cell type:", cell_type))
  }
  
  seu_subset <- subset(seu_obj, cells = cells)
  cat("Analyzing", length(cells), "cells for", cell_type, "\n")
  
  # Set identities to the condition we added earlier
  Idents(seu_subset) <- "condition"
  
  # Helper function to run FindMarkers and standardize column names
  # so that your downstream MA and Volcano plots work without changes
  get_contrast <- function(ident_1, ident_2) {
    res <- FindMarkers(
      seu_subset, 
      ident.1 = ident_1, 
      ident.2 = ident_2, 
      test.use = "wilcox",
      logfc.threshold = logfc_thresh, # Set to 0 to keep non-sig genes for plotting clouds
      min.pct = min_pct               # Filter genes expressed in <10% of cells
    )
    
    # Rename columns to match DESeq2 expectations
    res <- res %>%
      rownames_to_column(var = "gene") %>%
      dplyr::rename(log2FoldChange = avg_log2FC, padj = p_val_adj, pvalue = p_val) %>%
      column_to_rownames(var = "gene")
    
    return(res)
  }
  
  cat("Running Wilcoxon contrasts...\n")
  
  # Extract contrasts (ident_1 vs ident_2)
  res_water_vs_ctrl <- get_contrast("MCT-Water", "Control")
  res_blum_vs_ctrl  <- get_contrast("MCT-Blumeria", "Control")
  res_blum_vs_water <- get_contrast("MCT-Blumeria", "MCT-Water")
  res_water_vs_blum <- get_contrast("MCT-Water", "MCT-Blumeria")
  
  # Print summaries (count of sig DEGs based on arbitrary threshold for console preview)
  cat("\n=== Wilcoxon Results Summary for", cell_type, "(padj < 0.05) ===\n")
  cat("Water vs Control DEGs:", sum(res_water_vs_ctrl$padj < 0.05, na.rm = TRUE), "\n")
  cat("Blumeria vs Control DEGs:", sum(res_blum_vs_ctrl$padj < 0.05, na.rm = TRUE), "\n")
  cat("Blumeria vs Water DEGs:", sum(res_blum_vs_water$padj < 0.05, na.rm = TRUE), "\n")
  cat("Water vs Blumeria DEGs:", sum(res_water_vs_blum$padj < 0.05, na.rm = TRUE), "\n\n")
  
  # Save results if requested
  results_list <- list(
    water_vs_ctrl = res_water_vs_ctrl,
    blum_vs_ctrl  = res_blum_vs_ctrl,
    blum_vs_water = res_blum_vs_water,
    water_vs_blum = res_water_vs_blum
  )
  
  if (save_results) {
    filename <- paste0("DEG_Wilcox_", gsub("[^A-Za-z0-9]", "_", cell_type), ".rds")
    saveRDS(results_list, file = paste0("03-analysis_scratch/", filename))
    cat("Results saved to:", filename, "\n")
  }
  
  return(results_list)
}

# ---- Run Wilcoxon on all cell types ----

deg_results <- list()
cell_types <- levels(seu_obj)  
for (cell_type in cell_types) {
  cat("\n=== Processing", cell_type, "===\n")
  deg_results[[cell_type]] <- run_wilcox_deg(seu_obj, cell_type, save_results = FALSE)
}


# Function to create one panel of MA plots
create_ma_plot <- function(deg_results, cell_types, contrast_name, 
                           padj_thresh = 0.05, lfc_thresh = 0.5) {
  
  # Collect data for all cell types
  plot_data <- data.frame()
  up_counts <- integer(length(cell_types))
  down_counts <- integer(length(cell_types))
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    
    # Extract the appropriate contrast (logic remains unbroken!)
    if (contrast_name == "MCT-Water vs Control") {
      res <- deg_results[[cell_type]]$water_vs_ctrl
    } else if (contrast_name == "MCT-Water vs MCT-Blumeria") {
      res <- deg_results[[cell_type]]$water_vs_blum
    } else if (contrast_name == "MCT-Blumeria vs Control") {
      res <- deg_results[[cell_type]]$blum_vs_ctrl
    }
    
    # Convert to data frame
    df <- as.data.frame(res) %>%
      mutate(
        gene = rownames(.),
        cell_type = cell_type,
        cell_type_num = i,  # x-axis position
        significant = !is.na(padj) & padj < padj_thresh & abs(log2FoldChange) > lfc_thresh,
        direction = case_when(
          significant & log2FoldChange > 0 ~ "UP",
          significant & log2FoldChange < 0 ~ "DOWN",
          TRUE ~ "NS"
        )
      )
    
    plot_data <- rbind(plot_data, df)
    
    # Count significant genes
    up_counts[i] <- sum(df$direction == "UP", na.rm = TRUE)
    down_counts[i] <- sum(df$direction == "DOWN", na.rm = TRUE)
  }
  
  # Set up cell type order and colors
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types)
  
  # ---- NEW: Format title for italics ----
  display_title <- gsub("Blumeria", "<i>Blumeria</i>", contrast_name)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = cell_type_num, y = log2FoldChange)) +
    ylim(-5,5)+
    # Plot non-significant points first (gray)
    geom_jitter(data = filter(plot_data, direction == "NS"),
                width = 0.25, height = 0, size = 0.3, alpha = 0.2, color = "gray70") +
    # Plot significant points on top
    geom_jitter(data = filter(plot_data, direction == "UP"),
                width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#E31A1C") +
    geom_jitter(data = filter(plot_data, direction == "DOWN"),
                width = 0.25, height = 0, size = 0.4, alpha = 0.6, color = "#1F78B4") +
    # Reference line at 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    # Labels and theme
    scale_x_continuous(breaks = 1:length(cell_types),
                       labels = cell_types) +
    labs(title = display_title, # <-- Use the formatted title here
         x = NULL,
         y = "Log2FoldChange") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      # ---- NEW: Use element_markdown to parse the <i> tags ----
      plot.title = ggtext::element_markdown(hjust = 0.5, face = "bold", size = 16) 
    )
  
  # Add count annotations at top and bottom
  
  y_max <- 5
  y_min <- -5
  
  # UP counts (top, in red)
  for (i in seq_along(cell_types)) {
    p <- p + annotate("text", x = i, 
                      y = y_max*0.95,
                      label = up_counts[i], color = "#E31A1C", size = 3, fontface = "bold")
  }
  
  # DOWN counts (bottom, in blue)
  for (i in seq_along(cell_types)) {
    p <- p + annotate("text", x = i, 
                      y = y_min*0.95,
                      label = down_counts[i], color = "#1F78B4", size = 3, fontface = "bold")
  }
  
  p <- p + annotate("text", x = length(cell_types) + 0.8, y = y_max * 0.95,
                    label = "UP", color = "#E31A1C", size = 4, fontface = "bold") +
    annotate("text", x = length(cell_types) + 0.8, y = y_min * 0.95,
             label = "DOWN", color = "#1F78B4", size = 4, fontface = "bold")
  
  return(p)
}


p_water_blum <- create_ma_plot(deg_results, rev(cell_types),padj_thresh = 0.05, "MCT-Water vs MCT-Blumeria")

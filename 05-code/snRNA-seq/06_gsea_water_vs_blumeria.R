# ==============================================================================
# Script: 06_gsea_water_vs_blumeria.R
# Purpose: Perform Gene Set Enrichment Analysis (GSEA) on pseudobulk DESeq2 
#          results comparing MCT-Water vs MCT-Blumeria.
# Inputs:  03-analysis_scratch/DEG_Cardiomyocytes.rds
#          03-analysis_scratch/DEG_Macrophages.rds
# Outputs: 04-results/CM_Water_vs_Blum_GSEA.csv
#          04-results/MAC_Water_vs_Blum_GSEA.csv
#          04-results/Figure_Pathway_Analysis_GSEA_WaterVsBlum.png
# ==============================================================================

# ---- 1. Load Required Libraries ----------------------------------------------
library(clusterProfiler)  # Core GSEA and pathway enrichment functions
library(enrichplot)       # Visualization tools for enrichment results
library(org.Rn.eg.db)     # Rat genome annotation database
library(dplyr)            # Data manipulation
library(stringr)          # Text wrapping and string handling
library(ggplot2)          # Publication-quality plotting
library(patchwork)        # Plot composition and dynamic alignment
library(cowplot)          # Additional layout and alignment tools
library(msigdbr)

# ---- 2. Input Data Loading ---------------------------------------------------
# Define all target cell types and load their pseudobulk DESeq2 results
cell_types <- c("Cardiomyocyte", "Fibroblast", "Macrophage","Monocyte","Dendritic Cell",
                "Neutrophil","T Cell", "B Cell", "NK Cell", "Pericyte", "Vascular EC",
                "Endocardial EC", "Lymphatic EC", "Neuronal", "Epicardial EC")

deg_data <- lapply(cell_types, function(ct) {
  file_path <- sprintf("03-analysis_scratch/DEG_%s.rds", ct)
  if (file.exists(file_path)) {
    readRDS(file_path)
  } else {
    warning(sprintf("File %s not found. Skipping %s.", file_path, ct))
    NULL
  }
})
names(deg_data) <- cell_types

# Align cache location and fetch Rat Hallmark gene sets
Sys.setenv(R_USER_CACHE_DIR = tempdir())
hallmark_t2g <- msigdbr(species = "Rattus norvegicus", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene)

# ---- 3. Helper Functions -----------------------------------------------------

#' Extract and Rank Gene List for GSEA
#' 
#' @param res DESeq2 results object or data frame containing gene statistics
#' @return A sorted, named vector of log2FoldChange values indexed by ENTREZID
get_ranked_list <- function(res) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  # Convert gene SYMBOLs to ENTREZIDs required by clusterProfiler
  id_map <- bitr(
    df$gene, 
    fromType = "SYMBOL", 
    toType   = "ENTREZID", 
    OrgDb    = org.Rn.eg.db
  )
  
  # Merge mapped IDs back to DEG results and drop unmapped/NA values
  df <- merge(df, id_map, by.x = "gene", by.y = "SYMBOL")
  df <- df[!is.na(df$ENTREZID) & !is.na(df$log2FoldChange), ]
  
  # Construct named vector ranked by Log2 Fold Change
  gene_list <- df$log2FoldChange
  names(gene_list) <- df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

#' Run MSigDB Hallmark GSEA
#' 
#' @param ranked_vec Named vector of ranked gene statistics (ENTREZIDs)
#' @param t2g TERM2GENE dataframe from msigdbr
#' @return gseaResult object containing enriched Hallmark terms
run_hallmark_gsea <- function(ranked_vec, t2g) {
  set.seed(42) # Ensure reproducible random permutations
  
  res <- GSEA(
    geneList      = ranked_vec,
    TERM2GENE     = t2g,
    pvalueCutoff  = 0.05,            # FDR significance threshold
    pAdjustMethod = "BH",              # Benjamini-Hochberg adjustment
    minGSSize     = 10,
    maxGSSize     = 500,
    verbose       = FALSE
  )
  
  # Clean up formatting for plotting (e.g., "HALLMARK_GLYCOLYSIS" -> "Glycolysis")
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res@result$Description <- res@result$ID %>%
      gsub("HALLMARK_", "", .) %>%
      gsub("_", " ", .) %>%
      stringr::str_to_title()
  }
  
  return(res)
}

#' Plot Top Enriched GSEA Pathways as Horizontal Bar Charts
#' 
#' @param gsea_result GSEA result object (converted to gene symbols)
#' @param top_n Number of top pathways to display
#' @param direction Enriched direction ("up" for NES > 0, "down" for NES < 0)
#' @param drop_keywords Vector of regex terms to filter non-cardiac pathways
#' @param plot_title Title header for the subplot panel
#' @return A ggplot bar chart object
plot_gsea <- function(gsea_result, top_n = NULL, direction = "up", 
                      drop_keywords = NULL, plot_title = NULL) { 
  
  df <- as.data.frame(gsea_result)
  if (nrow(df) == 0) return(NULL)
  
  # Remove unwanted keywords if specified
  if (!is.null(drop_keywords)) {
    pattern <- paste(drop_keywords, collapse = "|")
    df <- df[!grepl(pattern, df$Description, ignore.case = TRUE), ]
  }
  
  # Filter by NES direction
  if (direction == "up") {
    df <- subset(df, NES > 0)
    bar_color <- "#FC8D62"
  } else if (direction == "down") {
    df <- subset(df, NES < 0)
    bar_color <- "#8DA0CB"
  }
  
  # Exit early if no pathways match the direction filter
  if (nrow(df) == 0) return(NULL)
  
  # Set dynamic axis limits after confirming data exists
  if (direction == "up") {
    x_limits <- c(0, max(3, ceiling(max(df$NES, na.rm = TRUE))))
  } else {
    x_limits <- c(min(-3, floor(min(df$NES, na.rm = TRUE))), 0)
  }
  
  # Select pathways
  plot_df <- df %>% arrange(p.adjust)
  if (!is.null(top_n)) {
    plot_df <- plot_df %>% slice_head(n = top_n)
  }
  
  # Order factors by NES so bars plot in order of magnitude
  plot_df <- plot_df %>%
    arrange(NES) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # Generate plot
  ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2, fill = bar_color) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) + 
    scale_x_continuous(limits = x_limits) +
    theme_bw(base_size = 14) + 
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 12, color = "black", lineheight = 0.8, face = "bold"),
      axis.text.x        = element_text(size = 12, color = "black"),
      axis.title.x       = element_text(size = 10, face = "bold", margin = margin(t = 5)),
      axis.title.y       = element_blank(),
      plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin        = margin(t = 15, r = 15, b = 15, l = 10)
    ) +
    labs(x = NULL, title = plot_title)
}

# ---- 4. Define Filtering Blacklist -------------------------------------------
# Keywords to strip non-cardiovascular neuronal/synaptic noise from GO annotation
# cardiac_blacklist <- c("behavior", "axon", "synapse", "neurotrans", "postsynaptic", 
                       # "ranvier", "AMPA", "glutamate", "dopamine", "presynaptic")

# ---- 5. Execute GSEA & Generate Figures Across Cell Types ---------------------
gsea_results_list <- list()
all_plots         <- list()

for (ct in names(deg_data)) {
  if (is.null(deg_data[[ct]])) next
  
  message(sprintf("Processing Hallmark GSEA for %s...", ct))
  
  # 1. Extract contrast and run Hallmark GSEA
  ranked_list <- get_ranked_list(deg_data[[ct]]$water_vs_blum)
  gsea_res    <- run_hallmark_gsea(ranked_list, hallmark_t2g)
  
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_res <- setReadable(gsea_res, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    gsea_results_list[[ct]] <- gsea_res
    
    # Clean filename by replacing spaces with underscores
    clean_ct_name <- gsub(" ", "_", toupper(ct))
    
    # 2. Export CSV Table
    write.csv(
      as.data.frame(gsea_res), 
      sprintf("04-results/%s_Water_vs_Blum_Hallmark_GSEA.csv", clean_ct_name), 
      row.names = FALSE
    )
    
    # Generate Activated and Suppressed Plots for ALL significant terms
    p_up   <- plot_gsea(gsea_res, top_n = NULL, direction = "up", 
                        plot_title = sprintf("%s: Activated Pathways", ct))
    
    p_down <- plot_gsea(gsea_res, top_n = NULL, direction = "down", 
                        plot_title = sprintf("%s: Suppressed Pathways", ct))
    
    # Safely collect non-empty plots
    if (!is.null(p_up) && !is.null(p_up$data) && nrow(p_up$data) > 0) {
      all_plots[[paste0(ct, "_up")]] <- p_up
    }
    if (!is.null(p_down) && !is.null(p_down$data) && nrow(p_down$data) > 0) {
      all_plots[[paste0(ct, "_down")]] <- p_down
    }
  } else {
    message(sprintf("No significant Hallmark pathways found for %s.", ct))
  }
}

# ---- 6. Compose and Save Final Master Figure ---------------------------------
if (length(all_plots) > 0) {
  message("Assembling final publication figure...")
  
  # Extract the number of pathways (rows) in each plot to use as relative heights
  panel_heights <- sapply(all_plots, function(p) nrow(p$data))
  
  # Stack all non-empty subplots vertically and apply proportional heights
  final_figure <- wrap_plots(all_plots, ncol = 1) +
    plot_layout(heights = panel_heights) +
    plot_annotation(
      title = "MCT-Water vs MCT-Blumeria",
      theme = theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5))
    )
  
  # Calculate a larger dynamic height (approx. 0.3 inches per pathway + margins)
  total_pathways <- sum(panel_heights)
  calc_height <- max(8, total_pathways * 0.3 + 3)
  
  # Save as PNG
  ggsave(
    filename = "04-results/Figure_Pathway_Analysis_GSEA_WaterVsBlum.png", 
    plot     = final_figure, 
    width    = 10, 
    height   = calc_height, 
    units    = "in", 
    dpi      = 300
  )
  
  # Save as PDF for publication
  ggsave(
    filename = "04-results/Figure_Pathway_Analysis_GSEA_WaterVsBlum.pdf", 
    plot     = final_figure, 
    width    = 10, 
    height   = calc_height, 
    units    = "in"
  )
}

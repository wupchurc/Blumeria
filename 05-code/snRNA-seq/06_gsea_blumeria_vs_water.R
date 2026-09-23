# ==============================================================================
# Script: 06_gsea_water_vs_blumeria.R
# Purpose: Perform Gene Set Enrichment Analysis (GSEA) on pseudobulk DESeq2 
#          results comparing MCT-Water vs MCT-Blumeria.
# Inputs:  03-analysis_scratch/DEG_Cardiomyocytes.rds
#          03-analysis_scratch/DEG_Macrophages.rds
# Outputs: 04-results/CM_Blum_vs_Water_GSEA.csv
#          04-results/MAC_Blum_vs_Water_GSEA.csv
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
cell_types <- c(
  "Cardiomyocyte", "Fibroblast", "Macrophage", "Monocyte",
  "Dendritic Cell", "Neutrophil", "T Cell", "B Cell", "NK Cell",
  "Pericyte", "Vascular EC", "Endocardial EC", "Lymphatic EC",
  "Neuronal"
)

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
#' Helper to wrap and left-pad y-axis labels to equalize plot panel widths
format_y_labels <- function(x, target_width = 36) {
  sapply(x, function(label) {
    wrapped <- stringr::str_wrap(label, width = target_width)
    lines   <- unlist(strsplit(wrapped, "\n"))
    padded  <- stringr::str_pad(lines, width = target_width, side = "left")
    paste(padded, collapse = "\n")
  })
}

#' Plot Top Enriched GSEA Pathways as Horizontal Bar Charts
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
  
  if (nrow(df) == 0) return(NULL)
  
  # Set dynamic axis limits
  if (direction == "up") {
    x_limits <- c(0, max(3, ceiling(max(df$NES, na.rm = TRUE))))
  } else {
    x_limits <- c(min(-3, floor(min(df$NES, na.rm = TRUE))), 0)
  }
  
  # Select and order pathways
  plot_df <- df %>% arrange(p.adjust)
  if (!is.null(top_n)) {
    plot_df <- plot_df %>% slice_head(n = top_n)
  }
  
  plot_df <- plot_df %>%
    arrange(NES) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # Generate plot
  ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_col(width = 0.5, color = "white", linewidth = 0.2, fill = bar_color) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    scale_y_discrete(
      labels = format_y_labels,
      expand = expansion(add = c(0.8, 0.8))
    ) + 
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_bw(base_size = 14) + 
    theme(
      plot.title.position = "plot", 
      panel.grid.major.y  = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.text.y         = element_text(size = 11, color = "black", lineheight = 0.8, face = "bold"),
      axis.text.x         = element_text(size = 10, color = "black"),
      axis.title.x        = element_text(size = 10, face = "plain", color = "grey20", margin = margin(t = 8)), # Cleaner, non-bold label
      axis.title.y        = element_blank(),
      plot.title          = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.margin         = margin(t = 15, r = 20, b = 15, l = 15)
    ) +
    labs(x = "Normalized Enrichment Score (NES)", title = plot_title)
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
  ranked_list <- get_ranked_list(deg_data[[ct]]$blum_vs_water)
  gsea_res    <- run_hallmark_gsea(ranked_list, hallmark_t2g)
  
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_res <- setReadable(gsea_res, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    gsea_results_list[[ct]] <- gsea_res
    
    clean_ct_name <- gsub(" ", "_", toupper(ct))
    
    # Export CSV Table
    write.csv(
      as.data.frame(gsea_res), 
      sprintf("04-results/%s_Blum_vs_Water_Hallmark_GSEA.csv", clean_ct_name), 
      row.names = FALSE
    )
    
    # Generate Plots
    p_up   <- plot_gsea(gsea_res, top_n = NULL, direction = "up", 
                        plot_title = sprintf("%s: Activated Pathways", ct))
    
    p_down <- plot_gsea(gsea_res, top_n = NULL, direction = "down", 
                        plot_title = sprintf("%s: Suppressed Pathways", ct))
    
    # Save UP plot
    if (!is.null(p_up) && !is.null(p_up$data) && nrow(p_up$data) > 0) {
      all_plots[[paste0(ct, "_up")]] <- p_up
      
      up_height <- (nrow(p_up$data) * 0.35) + 1.4
      
      ggsave(
        filename = sprintf("04-results/%s_Activated_Pathways.png", clean_ct_name),
        plot     = p_up,
        width    = 8.5,
        height   = up_height,
        units    = "in",
        dpi      = 300
      )
    }
    
    # Save DOWN plot
    if (!is.null(p_down) && !is.null(p_down$data) && nrow(p_down$data) > 0) {
      all_plots[[paste0(ct, "_down")]] <- p_down
      
      down_height <- (nrow(p_down$data) * 0.35) + 1.4
      
      ggsave(
        filename = sprintf("04-results/%s_Suppressed_Pathways.png", clean_ct_name),
        plot     = p_down,
        width    = 8.5,
        height   = down_height,
        units    = "in",
        dpi      = 300
      )
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
      title = "MCT-Blumeria vs MCT-Water",
      theme = theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5))
    )
  
  # Calculate a larger dynamic height (approx. 0.3 inches per pathway + margins)
  total_pathways <- sum(panel_heights)
  calc_height <- max(8, total_pathways * 0.3 + 3)
  
  # Save as PNG
  ggsave(
    filename = "04-results/Figure_Pathway_Analysis_GSEA_BlumVsWater.png", 
    plot     = final_figure, 
    width    = 10, 
    height   = calc_height, 
    units    = "in", 
    dpi      = 300
  )
  
  # Save as PDF for publication
  ggsave(
    filename = "04-results/Figure_Pathway_Analysis_GSEA_BlumVsWater.pdf", 
    plot     = final_figure, 
    width    = 10, 
    height   = calc_height, 
    units    = "in"
  )
}

# Summary Dotplot ----
# 1. Extract and combine pathways from GSEA results list
target_cts <- c("Cardiomyocyte", "Fibroblast", "Macrophage", "Pericyte")

desired_order <- c(
  # Inflammatory
  "Tnfa Signaling Via Nfkb", 
  "Interferon Gamma Response", 
  "Interferon Alpha Response", 
  "Inflammatory Response", 
  "Il6 Jak Stat3 Signaling", 
  "Il2 Stat5 Signaling",
  
  # Proliferation & Remodeling
  "G2m Checkpoint",
  "Epithelial Mesenchymal Transition",
  
  # Cell Stress & Survival
  "P53 Pathway",
  "Apoptosis",
  
  # Metabolism
  "Oxidative Phosphorylation",
  "Glycolysis"
)

df_combined <- bind_rows(lapply(target_cts, function(ct) {
  if (!is.null(gsea_results_list[[ct]])) {
    res <- as.data.frame(gsea_results_list[[ct]])
    res$Cell_Type <- ct
    return(res)
  }
})) %>%
  filter(Description %in% desired_order)

# 2. Map pathways to 4 distinct categories & set factor orders
df_faceted <- df_combined %>%
  mutate(Category = case_when(
    Description %in% c("Tnfa Signaling Via Nfkb", "Interferon Alpha Response", 
                       "Interferon Gamma Response", "Inflammatory Response",
                       "Il6 Jak Stat3 Signaling", "Il2 Stat5 Signaling") ~ "Inflammatory",
    
    Description %in% c("G2m Checkpoint", 
                       "Epithelial Mesenchymal Transition") ~ "Proliferation & Remodeling",
    
    Description %in% c("P53 Pathway", 
                       "Apoptosis") ~ "Cell Stress & Survival",
    
    Description %in% c("Oxidative Phosphorylation", 
                       "Glycolysis") ~ "Metabolism"
  )) %>%
  # Set Category factor levels (defines order of facet boxes from top to bottom)
  mutate(Category = factor(Category, levels = c(
    "Inflammatory", 
    "Proliferation & Remodeling", 
    "Cell Stress & Survival", 
    "Metabolism"
  ))) %>%
  # Set Description factor levels in reverse so level 1 is at the bottom of each facet
  mutate(Description = factor(Description, levels = rev(desired_order)))

# 3. Plot with 4 Facets
ggplot(df_faceted, aes(x = Cell_Type, y = Description, size = -log10(p.adjust), color = NES)) +
  geom_point() +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_color_gradient2(
    low = "#313695", 
    mid = "white", 
    high = "#A50026", 
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("-2.0", "-1.0", "0", "+1.0", "+2.0")
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    strip.text.y = element_text(face = "bold", angle = 0)
  ) +
  labs(x = NULL, y = NULL, title = NULL)

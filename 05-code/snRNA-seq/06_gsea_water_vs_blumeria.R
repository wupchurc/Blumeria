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

# ---- 2. Input Data Loading ---------------------------------------------------
# Load pseudobulk DESeq2 DEG list objects generated from script 05
cm_results  <- readRDS("03-analysis_scratch/DEG_Cardiomyocyte.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophage.rds")

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

#' Run Gene Ontology (Biological Process) GSEA
#' 
#' @param ranked_vec Named vector of ranked gene statistics (ENTREZIDs)
#' @return gseaResult object containing enriched GO terms
run_single_gsea <- function(ranked_vec) {
  set.seed(42) # Ensure reproducible random permutations
  
  gseGO(
    geneList      = ranked_vec,
    OrgDb         = org.Rn.eg.db,
    ont           = "BP",            # Target Biological Process ontology
    keyType       = "ENTREZID",
    pvalueCutoff  = 0.05,            # FDR significance threshold
    pAdjustMethod = "BH",              # Benjamini-Hochberg adjustment
    minGSSize     = 10,              # Exclude overly specific gene sets
    maxGSSize     = 500,             # Exclude broad, non-specific gene sets
    verbose       = FALSE,
    seed          = TRUE
  )
}

#' Plot Top Enriched GSEA Pathways as Horizontal Bar Charts
#' 
#' @param gsea_result GSEA result object (converted to gene symbols)
#' @param top_n Number of top pathways to display
#' @param direction Enriched direction ("up" for NES > 0, "down" for NES < 0)
#' @param drop_keywords Vector of regex terms to filter non-cardiac pathways
#' @param plot_title Title header for the subplot panel
#' @return A ggplot bar chart object
plot_gsea <- function(gsea_result, top_n = 10, direction = "up", 
                      drop_keywords = NULL, plot_title = NULL) { 
  
  df <- as.data.frame(gsea_result)
  
  # Remove irrelevant non-cardiac off-target pathways (e.g., neuronal terms)
  if (!is.null(drop_keywords)) {
    pattern <- paste(drop_keywords, collapse = "|")
    df <- df[!grepl(pattern, df$Description, ignore.case = TRUE), ]
  }
  
  # Filter by NES direction and assign color schemes
  if (direction == "up") {
    df <- subset(df, NES > 0)
    bar_color <- "#FC8D62" # Orange/Red for activated pathways
    x_limits  <- c(0, 3) 
  } else if (direction == "down") {
    df <- subset(df, NES < 0)
    bar_color <- "#8DA0CB" # Blue for suppressed pathways
    x_limits  <- c(-3, 0) 
  }
  
  # Select and order top pathways by adjusted p-value and NES magnitude
  plot_df <- df %>%
    arrange(p.adjust) %>%            
    slice_head(n = top_n) %>%        
    arrange(NES) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # Generate plot
  ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2, fill = bar_color) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) + # Wrap long pathway labels
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
cardiac_blacklist <- c("behavior", "axon", "synapse", "neurotrans", "postsynaptic", 
                       "ranvier", "AMPA", "glutamate", "dopamine", "presynaptic")

# ---- 5. Execute GSEA Workflow: Cardiomyocytes --------------------------------
message("Processing GSEA for Cardiomyocytes...")

# Extract contrast, calculate GSEA, and map IDs back to symbols
ranked_list_cm <- get_ranked_list(cm_results$water_vs_blum)
gsea_cm        <- run_single_gsea(ranked_list_cm)
gsea_cm        <- setReadable(gsea_cm, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")

# Generate Activated & Suppressed pathway plots
p_cm_up   <- plot_gsea(gsea_cm, top_n = 10, direction = "up", 
                       drop_keywords = cardiac_blacklist, 
                       plot_title = "Cardiomyocytes: Activated Pathways")

p_cm_down <- plot_gsea(gsea_cm, top_n = 10, direction = "down", 
                       drop_keywords = cardiac_blacklist, 
                       plot_title = "Cardiomyocytes: Suppressed Pathways")

# ---- 6. Execute GSEA Workflow: Macrophages -----------------------------------
message("Processing GSEA for Macrophages...")

# Extract contrast, calculate GSEA, and map IDs back to symbols
ranked_list_mac <- get_ranked_list(mac_results$water_vs_blum)
gsea_mac        <- run_single_gsea(ranked_list_mac)
gsea_mac        <- setReadable(gsea_mac, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")

# Generate Activated & Suppressed pathway plots
p_mac_up   <- plot_gsea(gsea_mac, top_n = 10, direction = "up", 
                        drop_keywords = cardiac_blacklist, 
                        plot_title = "Macrophages: Activated Pathways")

p_mac_down <- plot_gsea(gsea_mac, top_n = 10, direction = "down", 
                        drop_keywords = cardiac_blacklist, 
                        plot_title = "Macrophages: Suppressed Pathways")

# ---- 7. Export CSV Tabular Results -------------------------------------------
message("Exporting GSEA result tables...")

write.csv(as.data.frame(gsea_cm),  "04-results/CM_Water_vs_Blum_GSEA.csv",  row.names = FALSE)
write.csv(as.data.frame(gsea_mac), "04-results/MAC_Water_vs_Blum_GSEA.csv", row.names = FALSE)

# ---- 8. Compose and Save Final Figure ----------------------------------------
message("Assembling final publication figure...")

# Stack panels and dynamically scale row heights based on term count
final_figure <- p_cm_up / p_cm_down / p_mac_up / p_mac_down +
  plot_layout(heights = c(
    nrow(p_cm_up$data), 
    nrow(p_cm_down$data), 
    nrow(p_mac_up$data), 
    nrow(p_mac_down$data)
  ))

# Apply unified theme settings across all panels
final_figure <- final_figure &
  theme(
    plot.title = element_text(family = "Arial", size = 12, face = "bold", hjust = 0.5, margin = margin(b = 0)),
    plot.title.position = "panel"
  )

# Add master annotation title and shared X-axis label
final_figure <- final_figure +
  xlab("Normalized Enrichment Score (NES)") +
  plot_annotation(
    title = "MCT-Water vs MCT-Blumeria",
    theme = theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", hjust = 0.5))
  )

# Save figure file
ggsave(
  filename = "04-results/Figure_Pathway_Analysis_GSEA_WaterVsBlum.png", 
  plot     = final_figure, 
  width    = 11, 
  height   = 9, 
  units    = "in", 
  dpi      = 300
)

message("GSEA workflow completed successfully!")

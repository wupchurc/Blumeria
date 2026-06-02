# ---- Import Packages and Data ----
library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)
library(DESeq2)
library(ggtangle)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggtext)
library(patchwork)
library(cowplot)

cm_results <- readRDS("03-analysis_scratch/DEG_Cardiomyocytes.rds")
mac_results <- readRDS("03-analysis_scratch/DEG_Macrophages.rds")


# Function to get ranked list
get_ranked_list <- function(res) {
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  id_map <- bitr(df$gene, fromType="SYMBOL", 
                 toType="ENTREZID", OrgDb=org.Rn.eg.db)
  df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
  df <- df[!is.na(df$ENTREZID) & !is.na(df$log2FoldChange), ]
  gene_list <- df$log2FoldChange
  names(gene_list) <- df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}
# Set random seed
# set.seed(42) # Can be any number, 42 is a classic choice

run_single_gsea_all <- function(ranked_vec) {
  set.seed(42)
  gseGO(
    geneList = ranked_vec,
    OrgDb = org.Rn.eg.db,
    ont = "ALL",           # Changed from "BP" to "ALL"
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE,
    seed = TRUE
  )
}

plot_gsea <- function(gsea_result, top_n = 10, direction = "both", 
                      drop_keywords = NULL, target_ontology = "BP",
                      plot_title = NULL) { # Added a title parameter
  
  df <- as.data.frame(gsea_result)
  df <- subset(df, ONTOLOGY == target_ontology)
  
  if (!is.null(drop_keywords)) {
    pattern <- paste(drop_keywords, collapse = "|")
    df <- df[!grepl(pattern, df$Description, ignore.case = TRUE), ]
  }
  
  # Filter by direction
  if (direction == "up") {
    df <- subset(df, NES > 0)
    bar_color <- "#FC8D62" 
    x_limits <- c(0, 3) 
  } else if (direction == "down") {
    df <- subset(df, NES < 0)
    bar_color <- "#8DA0CB" 
    x_limits <- c(-3, 0) 
  }
  
  plot_df <- df %>%
    arrange(p.adjust) %>%            
    slice_head(n = top_n) %>%        
    arrange(NES) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2, fill = bar_color) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    
    # CHANGE 1: Use str_wrap instead of str_trunc to stack long pathway names cleanly
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    scale_x_continuous(limits = x_limits) +
    
    # CHANGE 2: Increase base_size for overall larger text
    theme_bw(base_size = 14) + 
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      # CHANGE 3: Adjust lineheight for wrapped text and increase font weight/size
      axis.text.y = element_text(size = 12, color = "black", lineheight = 0.8, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      
      # CHANGE 4: Add margins so long text strings aren't squeezed against the edge
      plot.margin = margin(t = 15, r = 15, b = 15, l = 10)
    ) +
    labs(x = "Normalized Enrichment Score (NES)", title = plot_title)
  
  return(p)
}

# ---- Cardiomyocytes ----

# 1. Extract the ranked list for your specific contrast
ranked_list_ctrl <- get_ranked_list(cm_results$water_vs_ctrl)
# 2. Run the GSEA across ALL ontologies
gsea_all_ctrl <- run_single_gsea_all(ranked_list_ctrl)
# 3. Convert ENTREZID back to SYMBOLs (so the core enrichment reads as gene names)
gsea_all_ctrl <- setReadable(gsea_all_ctrl, 
                             OrgDb = org.Rn.eg.db, 
                             keyType = "ENTREZID")

# 1. Extract the ranked list for your specific contrast
ranked_list_blum <- get_ranked_list(cm_results$water_vs_blum)
# 2. Run the GSEA across ALL ontologies
gsea_all_blum <- run_single_gsea_all(ranked_list_blum)
# 3. Convert ENTREZID back to SYMBOLs (so the core enrichment reads as gene names)
gsea_all_blum <- setReadable(gsea_all_blum, 
                             OrgDb = org.Rn.eg.db, 
                             keyType = "ENTREZID")

cardiac_blacklist <- c("behavior", "axon", "synapse", "neurotrans", "postsynaptic", 
                       "ranvier", "AMPA", "glutamate", "dopamine", "presynaptic")

# Water vs Blumeria
# 1. Generate the clean BP plots
p_wvb_suppressed <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "down", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

p_wvb_activated <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "up", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

# Water vs Blumeria
p_wvb_suppressed <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "down", 
  drop_keywords = cardiac_blacklist, target_ontology = "BP",
  plot_title = "Suppressed Pathways: Water vs Blumeria"
)

p_wvb_activated <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "up", 
  drop_keywords = cardiac_blacklist, target_ontology = "BP",
  plot_title = "Activated Pathways: Water vs Blumeria"
)

# 2. Count rows to automatically scale the panel heights flawlessly
rows_suppressed_c <- nrow(p_wvb_suppressed$data)
rows_activated_d  <- nrow(p_wvb_activated$data)

# Water vs Control
# 1. Generate the clean BP plots
p_wvc_suppressed <- plot_gsea(
  gsea_all_ctrl, top_n = 10, direction = "down", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)



p_wvc_activated <- plot_gsea(
  gsea_all_ctrl, top_n = 10, direction = "up", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

# 2. Count rows to automatically scale the panel heights flawlessly
rows_suppressed_a <- nrow(p_wvc_suppressed$data)
rows_activated_b  <- nrow(p_wvc_activated$data)

# 3. Stack vertically with proportional spacing
final_figure <- p_wvc_suppressed / p_wvc_activated / p_wvb_suppressed / p_wvb_activated + 
  plot_layout(heights = c(rows_suppressed_a, rows_activated_b,rows_suppressed_c, rows_activated_d )) 

final_figure <- final_figure +
  # plot_layout(heights = cm_heights, title = "Cardiomyocytes") +
  plot_annotation(
    title = "Cardiomyocytes", # Sub-header for the first stack
    tag_levels = list(c("A", "B", "C", "D"))
    ) &
  theme(
    plot.title = element_text(family = "Arial", size = 12, face = "bold", hjust = 0.5, margin = margin(b = -5)),
    plot.title.position = "panel",
    plot.tag = element_text(
      family = "Arial",   # Forces Arial font
      size = 17,          # Sets the precise font size
      face = "bold",      # Keeps it crisp and bold
      hjust = 0,          # Micro-adjustments to align perfectly top-left
      vjust = 1
    )
  )

final_figure_shrunk_cm <- final_figure & 
  theme(axis.text.y = element_text(size = 7, color = "black"))

ggsave(
  filename = "./04-results/Figure_Pathway_Analysis_GSEA_Cardiomyocytes.png",
  plot = final_figure_shrunk_cm,
  width = 11,      # 8 inches is great for a single-column layout
  height = 7,      # Shorter height works well since we removed CC and MF text
  units = "in",
  dpi = 300
)

gsea_list <- list(
  "Cardiomyocytes: Water vs Control" = gsea_all_ctrl,   # Replace with your actual object names
  "Cardiomyocytes: Water vs Blumeria" = gsea_all_blum
)
# Loop through each dataset and calculate the counts
for (dataset_name in names(gsea_list)) {
  gsea_df <- as.data.frame(gsea_list[[dataset_name]])
  
  # Filter for Biological Process (BP) and significant p.adjust
  sig_bp <- gsea_df %>% 
    filter(ONTOLOGY == "BP" & p.adjust < 0.05)
  
  # Split into Upregulated and Downregulated
  num_up <- sum(sig_bp$NES > 0)
  num_down <- sum(sig_bp$NES < 0)
  num_total <- nrow(sig_bp)
  
  # Print a publication-ready text summary to the console
  cat(paste0(
    "=== ", dataset_name, " ===\n",
    "Total significant BP pathways: ", num_total, "\n",
    "  - Activated (NES > 0): ", num_up, "\n",
    "  - Suppressed (NES < 0): ", num_down, "\n\n"
  ))
}

# ---- Macrophages ----
# 1. Extract the ranked list for your specific contrast
ranked_list_ctrl <- get_ranked_list(mac_results$water_vs_ctrl)
# 2. Run the GSEA across ALL ontologies
gsea_all_ctrl <- run_single_gsea_all(ranked_list_ctrl)
# 3. Convert ENTREZID back to SYMBOLs (so the core enrichment reads as gene names)
gsea_all_ctrl <- setReadable(gsea_all_ctrl, 
                             OrgDb = org.Rn.eg.db, 
                             keyType = "ENTREZID")

# 1. Extract the ranked list for your specific contrast
ranked_list_blum <- get_ranked_list(mac_results$water_vs_blum)
# 2. Run the GSEA across ALL ontologies
gsea_all_blum <- run_single_gsea_all(ranked_list_blum)
# 3. Convert ENTREZID back to SYMBOLs (so the core enrichment reads as gene names)
gsea_all_blum <- setReadable(gsea_all_blum, 
                             OrgDb = org.Rn.eg.db, 
                             keyType = "ENTREZID")

cardiac_blacklist <- c("behavior", "axon", "synapse", "neurotrans", "postsynaptic", 
                       "ranvier", "AMPA", "glutamate", "dopamine", "presynaptic")

# Water vs Blumeria
# 1. Generate the clean BP plots
p_wvb_suppressed <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "down", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

p_wvb_activated <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "up", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

# Water vs Blumeria
p_wvb_suppressed <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "down", 
  drop_keywords = cardiac_blacklist, target_ontology = "BP",
  plot_title = "Suppressed Pathways: Water vs Blumeria"
)

p_wvb_activated <- plot_gsea(
  gsea_all_blum, top_n = 10, direction = "up", 
  drop_keywords = cardiac_blacklist, target_ontology = "BP",
  plot_title = "Activated Pathways: Water vs Blumeria"
)

# 2. Count rows to automatically scale the panel heights flawlessly
rows_suppressed_c <- nrow(p_wvb_suppressed$data)
rows_activated_d  <- nrow(p_wvb_activated$data)

# Water vs Control
# 1. Generate the clean BP plots
p_wvc_suppressed <- plot_gsea(
  gsea_all_ctrl, top_n = 10, direction = "down", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

p_wvc_activated <- plot_gsea(
  gsea_all_ctrl, top_n = 10, direction = "up", drop_keywords = cardiac_blacklist, target_ontology = "BP"
)

# 2. Count rows to automatically scale the panel heights flawlessly
rows_suppressed_a <- nrow(p_wvc_suppressed$data)
rows_activated_b  <- nrow(p_wvc_activated$data)

# 3. Stack vertically with proportional spacing
final_figure <- p_wvc_suppressed / p_wvc_activated / p_wvb_suppressed / p_wvb_activated + 
  plot_layout(heights = c(rows_suppressed_a, rows_activated_b,rows_suppressed_c, rows_activated_d )) 

final_figure <- final_figure +
  plot_annotation(
    title = "Macrophages", # Sub-header for the second stack
    tag_levels = list(c("E", "F", "G", "H"))
    ) &
  theme(
    plot.title = element_text(family = "Arial", size = 12, face = "bold", hjust = 0.5, margin = margin(b = -5)),
    plot.title.position = "panel",
    plot.tag = element_text(
      family = "Arial",   # Forces Arial font
      size = 17,          # Sets the precise font size
      face = "bold",      # Keeps it crisp and bold
      hjust = 0,          # Micro-adjustments to align perfectly top-left
      vjust = 1
    )
  )

final_figure_shrunk_mac <- final_figure & 
  theme(axis.text.y = element_text(size = 7, color = "black"))

ggsave(
  filename = "./04-results/Figure_Pathway_Analysis_GSEA_Macrophages.png",
  plot = final_figure_shrunk_mac,
  width = 11,      # 8 inches is great for a single-column layout
  height = 7,      # Shorter height works well since we removed CC and MF text
  units = "in",
  dpi = 300
)

gsea_list <- list(
  "Macrophages: Water vs Control" = gsea_all_ctrl,   # Replace with your actual object names
  "Macrophages: Water vs Blumeria" = gsea_all_blum
)
# Loop through each dataset and calculate the counts
for (dataset_name in names(gsea_list)) {
  gsea_df <- as.data.frame(gsea_list[[dataset_name]])
  
  # Filter for Biological Process (BP) and significant p.adjust
  sig_bp <- gsea_df %>% 
    filter(ONTOLOGY == "BP" & p.adjust < 0.05)
  
  # Split into Upregulated and Downregulated
  num_up <- sum(sig_bp$NES > 0)
  num_down <- sum(sig_bp$NES < 0)
  num_total <- nrow(sig_bp)
  
  # Print a publication-ready text summary to the console
  cat(paste0(
    "=== ", dataset_name, " ===\n",
    "Total significant BP pathways: ", num_total, "\n",
    "  - Activated (NES > 0): ", num_up, "\n",
    "  - Suppressed (NES < 0): ", num_down, "\n\n"
  ))
}





# ----

# 1. Apply your x-limits and text shrinkage to both figures
# This guarantees they both span exactly from -3 to 3 on the x-axis
final_cardiomyocytes_prepped <- final_figure_shrunk_cm & 
  scale_x_continuous(limits = c(-3, 3)) &
  theme(axis.text.y = element_text(size = 7, color = "black"))

final_macrophages_prepped <- final_figure_shrunk_mac & 
  scale_x_continuous(limits = c(-3, 3)) &
  theme(axis.text.y = element_text(size = 7, color = "black"))

aligned_list <- cowplot::align_plots(final_cardiomyocytes_prepped, final_macrophages_prepped, align = "h", axis = "lr")

# 3. Save the Cardiomyocytes Figure
ggsave(
  filename = "./04-results/Figure_Pathway_Analysis_GSEA_Cardiomyocytes.png",
  plot = aligned_plots[[1]], # Grab the first aligned figure
  width = 11,      
  height = 7,      
  units = "in",
  dpi = 300
)

# 4. Save the Macrophages Figure
ggsave(
  filename = "./04-results/Figure_Pathway_Analysis_GSEA_Macrophages.png",
  plot = aligned_plots[[2]], # Grab the second aligned figure
  width = 11,      
  height = 7,      
  units = "in",
  dpi = 300
)




# ---- Combined GSEA ----
analyze_celltype_gsea <- function(results_obj, cell_name) {
  
  # Function to get ranked list
  get_ranked_list <- function(res) {
    df <- as.data.frame(res)
    df$gene <- rownames(res)
    id_map <- bitr(df$gene, fromType="SYMBOL", 
                   toType="ENTREZID", OrgDb=org.Rn.eg.db)
    df <- merge(df, id_map, by.x="gene", by.y="SYMBOL")
    df <- df[!is.na(df$ENTREZID) & !is.na(df$log2FoldChange), ]
    gene_list <- df$log2FoldChange
    names(gene_list) <- df$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    return(gene_list)
  }
  
  ranked_lists <- list(
    water_vs_ctrl = get_ranked_list(results_obj$water_vs_ctrl),
    water_vs_blum = get_ranked_list(results_obj$water_vs_blum)
  )
  
  gsea_df <- do.call(rbind, lapply(names(ranked_lists), function(nm) {
    data.frame(
      ENTREZID = names(ranked_lists[[nm]]),
      log2FoldChange = ranked_lists[[nm]],
      contrast_name = nm,
      stringsAsFactors = FALSE
    )
  }))
  # Set desired order here
  gsea_df$contrast_name <- factor(gsea_df$contrast_name,
                                  levels = c("water_vs_ctrl", "water_vs_blum"))
  
  comp_gsea <- compareCluster(
    ENTREZID | log2FoldChange ~ contrast_name,
    data = gsea_df,
    fun = "gseGO",
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
  
  # Make sure Cluster is ordered consistently for plotting
  comp_gsea@compareClusterResult$Cluster <- factor(
    comp_gsea@compareClusterResult$Cluster,
    levels = c("water_vs_ctrl", "water_vs_blum")
  )
  
  p_gsea <- dotplot(comp_gsea, showCategory = 5, font.size = 8,
                    split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste(cell_name, "- GO (BP)"))
  
  print(p_gsea)
  
  list(gsea = comp_gsea, ranked_lists = ranked_lists)
}

# ---- Single GSEA ----
run_single_gsea <- function(ranked_vec) {
  gseGO(
    geneList = ranked_vec,
    OrgDb = org.Rn.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}

make_single_gsea_list <- function(gsea_obj) {
  list(
    water_vs_ctrl = run_single_gsea(gsea_obj$ranked_lists$water_vs_ctrl),
    water_vs_blum = run_single_gsea(gsea_obj$ranked_lists$water_vs_blum)
  )
}

# ---- Cardiomyocytes ----
gsea_res <- analyze_celltype_gsea(cm_results, "Cardiomyocytes")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)
write.csv(
  gsea_df,
  file = "04-results/CM_GSEA_GO.csv",
  row.names = FALSE
)


# ---- Macrophages ----

gsea_res <- analyze_celltype_gsea(mac_results, "Macrophages")

gsea_df <- as.data.frame(gsea_res$gsea)
head(gsea_df)
write.csv(
  gsea_df,
  file = "04-results/MAC_GSEA_GO.csv",
  row.names = FALSE
)



# ---- Dotplots ----
# Simplify the combined result
gsea_simplified <- simplify(gsea_res$gsea, cutoff = 0.7, by = "p.adjust", 
                            measure = "Wang")
gsea_simplified_df <- as.data.frame(gsea_simplified)
write.csv(
  gsea_simplified_df,
  file = "04-results/CM_GSEA_compareCluster_results_simplified.csv",
  row.names = FALSE
)

# simplified dotplot
dotplot(gsea_simplified, showCategory = 5, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Macrophages") +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10)
  ) +
  labs(x = "", y = "")



# two separate plots for activated and suppressed
# Create a copy for activated only
gsea_activated <- gsea_simplified
gsea_activated@compareClusterResult <- gsea_activated@compareClusterResult[
  gsea_activated@compareClusterResult$NES > 0, ]

# Create a copy for suppressed only
gsea_suppressed <- gsea_simplified
gsea_suppressed@compareClusterResult <- gsea_suppressed@compareClusterResult[
  gsea_suppressed@compareClusterResult$NES < 0, ]

# Plot each
p1 <- dotplot(gsea_activated, showCategory = 5) +
  ggtitle("Activated Pathways") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10)) +
  labs(x = '')

p2 <- dotplot(gsea_suppressed, showCategory = 5) +
  ggtitle("Suppressed Pathways") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10)) +
  labs(x = '')

# Stack vertically
library(patchwork)
p1 / p2





# non simplified dot plot
dotplot(gsea_res$gsea, showCategory = 15, split = ".sign", size = "Count") + 
  facet_grid(. ~ .sign) + 
  ggtitle("CM - GSEA") +
  theme(axis.text.y = element_text(size = 8))

# ---- single plots ---- 
gsea_single <- make_single_gsea_list(gsea_res)
# Convert ENTREZID → SYMBOL
gsea_single$water_vs_ctrl <- setReadable(gsea_single$water_vs_ctrl, 
                                         OrgDb = org.Rn.eg.db, 
                                         keyType = "ENTREZID")
gsea_single$water_vs_blum <- setReadable(gsea_single$water_vs_blum, 
                                         OrgDb = org.Rn.eg.db, 
                                         keyType = "ENTREZID")

dotplot(gsea_single$water_vs_ctrl, showCategory = 10, split = ".sign") + ggtitle("CM water_vs_ctrl")
dotplot(gsea_single$water_vs_blum, showCategory = 10) + ggtitle("CM water_vs_blum")








# ---- Single GSEA ALL Ontologies----

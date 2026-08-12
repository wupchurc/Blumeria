# library(tidyverse)
library(microshades)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr) 
library(ggtext)


# Function ----
make_stacked_plot <- function(raw, sample_map_df, out_png, out_csv, my_palette, top_n = 10) {
  names(raw)[1] <- "Taxon"
  sample_cols <- names(raw)[-1]
  
  sample_info <- tibble(Sample = sample_cols) %>%
    mutate(Group = case_when(
      str_detect(Sample, "Control") ~ "Control",
      str_detect(Sample, "MCT_Water") ~ "MCT-Water",
      # str_detect(Sample, "MCT_Blumeria") ~ "MCT-Blumeria",
      str_detect(Sample, "MCT_Blumeria") ~ "MCT-<i>Blumeria</i>",
      TRUE ~ "Other"
    ))
  
  raw[sample_cols] <- lapply(raw[sample_cols], as.numeric)
  
  long <- raw %>%
    pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sample_info, by = "Sample") %>%
    left_join(sample_map_df, by = "Sample")
  
  top_taxa <- long %>%
    group_by(Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon)
  
  plot_df <- long %>%
    mutate(Taxon = if_else(Taxon %in% top_taxa, Taxon, "Other")) %>%
    group_by(Group, Sample, ShortName, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # plot_df$Group <- factor(plot_df$Group, levels = c("Control", "MCT-Water", "MCT-Blumeria"))
  # plot_df$Taxon <- factor(plot_df$Taxon, levels = c(rev(top_taxa), "Other"))
  # --- UPDATE: Updated factor levels to match markdown ---
  plot_df$Group <- factor(plot_df$Group, levels = c("Control", "MCT-Water", "MCT-<i>Blumeria</i>"))
  plot_df$Taxon <- factor(plot_df$Taxon, levels = c(rev(top_taxa), "Other"))
  
  # group_order <- c("Control", "MCT-Water", "MCT-Blumeria")
  group_order <- c("Control", "MCT-Water", "MCT-<i>Blumeria</i>")
  sample_order <- sample_map_df %>%
    mutate(Group = case_when(
      str_detect(Sample, "^Control") ~ "Control",
      str_detect(Sample, "^MCT_Water") ~ "MCT-Water",
      # str_detect(Sample, "^MCT_Blumeria") ~ "MCT-Blumeria",
      str_detect(Sample, "^MCT_Blumeria") ~ "MCT-<i>Blumeria</i>",
      TRUE ~ "Other"
    )) %>%
    arrange(factor(Group, levels = group_order), ShortName) %>%
    pull(ShortName)
  
  plot_df$ShortName <- factor(plot_df$ShortName, levels = sample_order)
  
  p <- ggplot(plot_df, aes(x = ShortName, y = Abundance, fill = Taxon)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    labs(
      title = "",
      x = NULL,
      y = "Relative abundance (%)",
      fill = "Genus"
    ) +
    scale_fill_manual(values = my_palette) +
    theme_bw(base_size = 12) +
    # guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey60"),
      # strip.text = element_text(size = 16, face = "bold"),
      # --- UPDATE: Change element_text to element_markdown ---
      strip.text = element_markdown(size = 16, face = "bold"),
      
      # --- CHANGES MADE HERE ---
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(), 
      # -------------------------
      
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.box = "horizontal",
      # legend.key.height = unit(0.4, "cm"),
      # legend.key.width = unit(0.7, "cm"),
      # legend.spacing.y = unit(0.1, "cm"),
      # legend.text = element_text(size = 9),
      # legend.title = element_text(size = 10),
      # plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(1.0, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  
  ggsave(out_png, p, width = 12, height = 7, dpi = 300, device = "png")
  write_csv(plot_df, out_csv)
  
  return(p)
}

# Metagenomics Call ----
meta_raw <- read.table(
  "04-results/metagenomics_merged_abundance_table_genus.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

meta_map <- tibble(
  Sample = c(
    "Control_SP_11-Control_SP_11-GCCATAGCCATA",
    "Control_SP_12-Control_SP_12-TCTGGTTCTGGT",
    "Control_SP_13-Control_SP_13-TGGCTTTGGCTT",
    "Control_SP_14-Control_SP_14-AACACTAACACT",
    "Control_SP_15-Control_SP_15-TCAGGATCAGGA",
    "Control_SP_16-Control_SP_16-TGGTCCTGGTCC",
    "Control_SP_17-Control_SP_17-CAACTGCAACTG",
    "MCT_Blumeria_1-MCT_Blumeria_1-AACACTAACACT",
    "MCT_Blumeria_3-MCT_Blumeria_3-TCAGGATCAGGA",
    "MCT_Blumeria_4-MCT_Blumeria_4-TGGTCCTGGTCC",
    "MCT_Blumeria_6-MCT_Blumeria_6-CAACTGCAACTG",
    "MCT_Blumeria_7-MCT_Blumeria_7-CGGACCCGGACC",
    "MCT_Blumeria_8-MCT_Blumeria_8-ATCGAGATCGAG",
    "MCT_Blumeria_9-MCT_Blumeria_9-ATGGTGATGGTG",
    "MCT_Blumeria_10-MCT_Blumeria_10-CTTAAGCTTAAG",
    "MCT_Blumeria_11-MCT_Blumeria_11-GAGTGCGAGTGC",
    "MCT_Water_SP_1-MCT_Water_SP_1-CGGACCCGGACC",
    "MCT_Water_SP_2-MCT_Water_SP_2-ATCGAGATCGAG",
    "MCT_Water_SP_5-MCT_Water_SP_5-ATGGTGATGGTG",
    "MCT_Water_SP_6-MCT_Water_SP_6-CTTAAGCTTAAG",
    "MCT_Water_SP_9-MCT_Water_SP_9-GAGTGCGAGTGC",
    "MCT_Water_SP_10-MCT_Water_SP_10-GCCATAGCCATA",
    "MCT_Water_SP_11-MCT_Water_SP_11-TCTGGTTCTGGT",
    "MCT_Water_SP_13-MCT_Water_SP_13-TGGCTTTGGCTT"
  ),
  ShortName = c(
    "C11","C12","C13","C14","C15","C16","C17",
    "B1","B3","B4","B6","B7","B8","B9","B10","B11",
    "W1","W2","W5","W6","W9","W10","W11","W13"
  )
)

meta_palette <- c(
  "Lactobacillus" = "#D7B5D8",
  "Ligilactobacillus" = "#B8A1E3",
  "Xylanibacter" = "#9CC9F5",
  "Muribaculum" = "#7EC8D8",
  "GGB14010" = "#6AB7A8",
  "Lachnospiraceae_unclassified" = "#8DCB9D",
  "GGB1515" = "#B9D98A",
  "Akkermansia" = "#D9E27A",
  "Limosilactobacillus" = "#F3C88E",
  "Parabacteroides" = "#E7A6A6",
  "Other" = "#E6E6E6"
)

p_meta <- make_stacked_plot(
  raw = meta_raw,
  sample_map_df = meta_map,
  out_png = "04-results/stacked_bar_metagenomics.png",
  out_csv = "04-results/stacked_bar_metagenomics_fixed_data.csv",
  my_palette = meta_palette,
  top_n = 10
)

ggsave("stacked_bar_metagenomics.png", p_meta, width = 12, height = 7, dpi = 300, device="png")

# ITS2 Call ----
# its_raw <- read_csv("04-results/ITS2_taxa_genus_relative_abundance_fungi_only.txt", show_col_types = FALSE)

its_raw <- read_tsv(
  "04-results/ITS2_taxa_genus_relative_abundance_fungi_only.txt",
  show_col_types = FALSE,
  name_repair = "unique"
)

its_raw <- its_raw %>%
  select(Genus, everything()) %>%
  rename(Taxon = Genus)

drop_cols <- c(2, (ncol(its_raw) - 14):ncol(its_raw))
its_raw <- its_raw[, -drop_cols]
its_raw$Taxon <- sub("^g__", "", its_raw$Taxon)
# names(its_raw)[1] <- "Taxon"
sample_cols <- names(its_raw)[-1]

its_raw[sample_cols] <- lapply(its_raw[sample_cols], function(x) as.numeric(x) * 100)


its_map <- tibble(
  Sample = c(
    "Control_SP_11_S17_",
    "Control_SP_12_S18_",
    "Control_SP_13_S19_",
    "Control_SP_14_S20_",
    "Control_SP_15_S21_",
    "Control_SP_16_S22_",
    "Control_SP_17_S23_",
    "MCT_Blumeria_1_S34_",
    "MCT_Blumeria_10_S32_",
    "MCT_Blumeria_11_S33_",
    "MCT_Blumeria_3_S35_",
    "MCT_Blumeria_4_S36_",
    "MCT_Blumeria_6_S37_",
    "MCT_Blumeria_7_S38_",
    "MCT_Blumeria_8_S39_",
    "MCT_Blumeria_9_S40_",
    "MCT_Water_SP_1_S44_",
    "MCT_Water_SP_10_S41_",
    "MCT_Water_SP_11_S42_",
    "MCT_Water_SP_13_S43_",
    "MCT_Water_SP_2_S45_",
    "MCT_Water_SP_5_S46_",
    "MCT_Water_SP_6_S47_",
    "MCT_Water_SP_9_S48_"
  ),
  ShortName = c(
    "C11","C12","C13","C14","C15","C16","C17",
    "B1","B10","B11","B3","B4","B6","B7","B8","B9",
    "W1","W10","W11","W13","W2","W5","W6","W9"
  )
)

its_palette <- c(
  "Thermomyces" = "#D7B5D8",
  "Myxotrichum" = "#B8A1E3",
  "Vishniacozyma" = "#9CC9F5",
  "Calophoma" = "#7EC8D8",
  "Oidiodendron" = "#6AB7A8",
  "Unassigned" = "#8DCB9D",
  "Talaromyces" = "#B9D98A",
  "Penicillium" = "#D9E27A",
  "Cladosporium" = "#F3C88E",
  "Aspergillus" = "#E7A6A6",
  "Other" = "#E6E6E6"
)

p_its <- make_stacked_plot(
  raw = its_raw,
  sample_map_df = its_map,
  out_png = "04-results/stacked_bar_ITS2.png",
  out_csv = "04-results/stacked_bar_ITS2_fungi_data.csv",
  my_palette = its_palette,
  top_n = 5
)

ggsave("stacked_bar_its2.png", p_its, width = 12, height = 7, dpi = 300, device="png")


library(tidyverse)

# 1. Load the FULL abundance table and the GTDB dictionary
abundance_full <- read_tsv("04-results/metagenomics_merged_abundance_table.txt", comment = "#")
gtdb_map <- read_tsv("/scratch.global/wupchurc/metaphlan_db/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv", 
                     col_names = c("sgb_id", "gtdb_taxonomy"))

# 2. Process the table, map the taxonomy, and add the clean genus column
gtdb_genus_table <- abundance_full %>%
  # Keep only the highly specific Species (SGB) rows
  filter(str_detect(clade_name, "s__.*SGB")) %>%
  
  # Extract the SGB ID to match the dictionary
  mutate(clean_sgb = str_extract(clade_name, "SGB[0-9]+")) %>%
  left_join(gtdb_map, by = c("clean_sgb" = "sgb_id")) %>%
  
  # Extract the GTDB lineage up to the Genus level
  mutate(gtdb_genus = str_extract(gtdb_taxonomy, "^[^;]+;[^;]+;[^;]+;[^;]+;[^;]+;g__[^;]+")) %>%
  mutate(gtdb_genus = str_replace_all(gtdb_genus, ";", "|")) %>%
  
  # Fall back to original genus if no GTDB map exists
  mutate(original_genus = str_extract(clade_name, "^.*\\|g__[^\\|]+")) %>%
  
  # Set the full taxonomic path as the 'clade_name'
  mutate(clade_name = coalesce(gtdb_genus, original_genus)) %>%
  
  # ---> NEW: Extract just the clean genus name without the "g__"
  mutate(clean_genus = str_extract(clade_name, "g__.*$")) %>%
  mutate(clean_genus = str_remove(clean_genus, "g__")) %>%
  
  # Drop the intermediate mapping columns
  select(-any_of("NCBI_tax_id"), -clean_sgb, -gtdb_taxonomy, -gtdb_genus, -original_genus) %>%
  
  # Reorder so 'clean_genus' is column #2, right next to 'clade_name'
  select(clade_name, clean_genus, everything()) %>%
  
  # Group by BOTH name columns so we don't lose the new column when summing the abundances
  group_by(clade_name, clean_genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

# 3. Save your finalized GTDB Genus table
# write_tsv(gtdb_genus_table, "04-results/metagenomics_merged_abundance_table_gtdb_genus.txt")
# 3. Save your finalized GTDB Genus table as a CSV
write_csv(gtdb_genus_table, "04-results/metagenomics_merged_abundance_table_gtdb_genus.csv")

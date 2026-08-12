library(readr)
library(dplyr)

df <- read_tsv("04-results/metagenomics_merged_abundance_table.txt",
               skip = 1)

genus_df <- df %>%
  filter(grepl("\\|g__", .[[1]]) & !grepl("\\|s__", .[[1]]) & !grepl("\\|t__", .[[1]]))
genus_df[[1]] <- sub("^.*\\|g__", "", genus_df[[1]])


# genus_df <- df %>%
  # filter(grepl("\\|g__", .[[1]]) & !grepl("\\|s__", .[[1]]) & !grepl("\\|t__", .[[1]])) %>%
  # mutate(genus = sub("^.*\\|g__", "", .[[1]])) %>%
  # select(genus, everything(), -all_of(names(df)[1])) %>%
  # group_by(genus) %>%
  # summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

write_tsv(genus_df, "04-results/metagenomics_merged_abundance_table_genus.txt")

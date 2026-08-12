library(tidyverse)

# 1. Load the data
data <- read_tsv("04-results/metagenomics_merged_abundance_table_genus.txt", show_col_types = FALSE)

# 2. Filter out "GGB" taxa and prep for row insertion
data_filtered <- data %>%
  # Filter out rows where clade_name starts with "GGB"
  filter(!str_starts(clade_name, "GGB")) %>%
  # Convert all columns to character so we can safely insert text labels later
  mutate(across(everything(), as.character))

# 3. Create the condition row based on column headers using case_when
col_names <- names(data_filtered)

condition_values <- case_when(
  col_names == "clade_name" ~ "Condition",
  str_starts(col_names, "Control") ~ "Control",
  str_starts(col_names, "MCT_Blumeria") ~ "MCT-Blumeria",
  str_starts(col_names, "MCT_Water") ~ "MCT-Water",
  TRUE ~ "Unknown"
)

# Convert the named vector into a 1-row tibble
condition_row <- bind_rows(setNames(condition_values, col_names))

# OPTION A: Insert condition row right below the header
final_data_a <- bind_rows(condition_row, data_filtered)

write_csv(final_data_a, "04-results/metagenomics_no_ggb_for_MetaboAnalyst.csv")

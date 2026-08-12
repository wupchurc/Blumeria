library(readr)
library(dplyr)
library(stringr)

its <- readr::read_tsv("ITS_combined_sequences_taxa.txt")

# Taxonomy columns are the last 7
tax_cols <- tail(names(its), 7)

# Sample columns = all columns except first (ID) and last 7 (taxonomy)
sample_cols <- names(its)[2:(ncol(its) - length(tax_cols))]

sample_cols

# Include every sequenced row ----

its <- its %>%
  mutate(
    Genus_clean = Genus,
    Genus_clean = if_else(is.na(Genus_clean) | Genus_clean == "NA",
                          "Unassigned", Genus_clean),
    Genus_clean = sub("^g__", "", Genus_clean)
  )

genus_counts <- its %>%
  group_by(Genus_clean) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop")

# Quick check
head(genus_counts)

# Compute per-sample relative abundance
genus_rel <- genus_counts %>%
  mutate(across(
    all_of(sample_cols),
    ~ .x / sum(.x),
    .names = "{.col}"          # keep same names
  ))

# Optional sanity check: columns should each sum to ~1
colSums(genus_rel[, sample_cols])

# Rename Genus_clean to Genus, and ensure it is the first column
output <- genus_rel %>%
  rename(Genus = Genus_clean) %>%
  select(Genus, all_of(sample_cols))

head(output)

# Build condition vector matched to sample names
condition_vec <- case_when(
  grepl("^Control",      sample_cols) ~ "Control",
  grepl("^MCT_Blumeria", sample_cols) ~ "MCT-Blumeria",
  grepl("^MCT_Water",    sample_cols) ~ "MCT-Water",
  TRUE ~ "Unknown"  # fallback if anything doesn't match
)

condition_row <- tibble(
  Genus = "Condition",
  !!!setNames(as.list(condition_vec), sample_cols)
)

# For export only: coerce sample columns to character
output_chr <- output %>%
  mutate(across(all_of(sample_cols), as.character))

# Bind condition row *above* all genera
export_mat <- bind_rows(condition_row, output_chr)

write_csv(export_mat, "ITS_genus_relative_abundance.csv")

# include only fungi sequences ----

its_fungi <- its %>%
  filter(str_detect(Kingdom, "Fungi"))

its_fungi <- its_fungi %>%
  mutate(
    Genus_clean = Genus,
    Genus_clean = if_else(is.na(Genus_clean) | Genus_clean == "NA",
                          "Unassigned", Genus_clean),
    Genus_clean = sub("^g__", "", Genus_clean)
  )

genus_counts <- its_fungi %>%
  group_by(Genus_clean) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop")

# Quick check
head(genus_counts)

# Compute per-sample relative abundance
genus_rel <- genus_counts %>%
  mutate(across(
    all_of(sample_cols),
    ~ .x / sum(.x),
    .names = "{.col}"          # keep same names
  ))

# Optional sanity check: columns should each sum to ~1
colSums(genus_rel[, sample_cols])

# Rename Genus_clean to Genus, and ensure it is the first column
output <- genus_rel %>%
  rename(Genus = Genus_clean) %>%
  select(Genus, all_of(sample_cols))

head(output)

# Build condition vector matched to sample names
condition_vec <- case_when(
  grepl("^Control",      sample_cols) ~ "Control",
  grepl("^MCT_Blumeria", sample_cols) ~ "MCT-Blumeria",
  grepl("^MCT_Water",    sample_cols) ~ "MCT-Water",
  TRUE ~ "Unknown"  # fallback if anything doesn't match
)

condition_row <- tibble(
  Genus = "Condition",
  !!!setNames(as.list(condition_vec), sample_cols)
)

# For export only: coerce sample columns to character
output_chr <- output %>%
  mutate(across(all_of(sample_cols), as.character))

# Bind condition row *above* all genera
export_mat <- bind_rows(condition_row, output_chr)

write_csv(export_mat, "ITS_fungi_genus_relative_abundance.csv")

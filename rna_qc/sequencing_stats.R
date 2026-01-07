# Load libraries ----------------------------------------------------------

library(scales)
library(tidyverse)

# Define variables --------------------------------------------------------


# Load data ---------------------------------------------------------------

# Find all files
files = list.files(
  path = "/Volumes/bcn_synaptopathy_in_ftd_als/synaptosome_long-read/data/nanostat",
  pattern = "barcode.*_nanostat_out\\.txt$",
  recursive = TRUE, # Search all subdirectories
  full.names = TRUE # Return full paths
)

# Create a function to read one file
read_sequencing_stats = function(file) {
  # Extract barcode name from file path
  sample = str_extract(file, "barcode[0-9]+")
  
  # Read tsv
  read_tsv(
    file,
    comment = "#",
    col_names = TRUE,
    show_col_types = FALSE # Suppress noise
  ) %>%
    mutate(
      sample = sample,
    ) %>%
    rename(
      stat = Metrics,
      value = dataset
    )
}

# Read all files and combine
sequencing_stats = files %>%
  # map() applies a function once per element and returns a list
  map(read_sequencing_stats) %>%
  # list_rbind() takes a list of dataframes/tibbles and stacks them row-wise
  list_rbind()

stat_order = c(
  "number_of_reads",
  "median_read_length",
  "n50",
  "longest_read_(with_Q):1",
  "median_qual"
)

# Pivot wider and clean
sequencing_stats_wide = sequencing_stats %>%
  pivot_wider(
    names_from = "sample",
    values_from = "value"
  ) %>%
  rename(
    WT_Cyt_1 = barcode01,
    Q331K_Cyt_1 = barcode02, 
    WT_Syn_1 = barcode03,
    Q331K_Syn_1 = barcode04,
    WT_Cyt_2 = barcode05,
    Q331K_Cyt_2 = barcode06,
    WT_Syn_2 = barcode07,
    Q331K_Syn_2 = barcode08,
    WT_Cyt_3 = barcode09,
    Q331K_Cyt_3 = barcode10,
    WT_Syn_3 = barcode11,
    Q331K_Syn_3 = barcode12,
    WT_Cyt_4 = barcode13,
    Q331K_Cyt_4 = barcode14,
    WT_Syn_4 = barcode15,
    Q331K_Syn_4 = barcode16
  ) %>%
  filter(stat %in% stat_order) %>%
  mutate(stat = factor(stat, levels=stat_order)) %>%
  arrange(stat) %>%
  select(stat,
         starts_with("WT_Cyt"),
         starts_with("Q331K_Cyt"),
         starts_with("WT_Syn"),
         starts_with("Q331K_Syn"))

# Reformat from scientific notation
sequencing_stats_formatted = sequencing_stats_wide %>%
  mutate(
    across(-stat, ~ case_when(
      stat == "longest_read_(with_Q):1" ~ as.character(.x),
      stat == "median_qual" ~ number(as.numeric(.x), accuracy=0.1),
      TRUE ~ comma(as.numeric(.x))
    ))
  )

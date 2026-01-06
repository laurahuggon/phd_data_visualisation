# Load libraries ----------------------------------------------------------

library(scales)
library(tidyverse)

# Define variables --------------------------------------------------------


# Load data ---------------------------------------------------------------

# Find all files
files = list.files(
  path = "/Volumes/bcn_synaptopathy_in_ftd_als/synaptosome_long-read/data/processed_no-annotation",
  pattern = "barcode.*_alignment_stats\\.tsv$",
  recursive = TRUE, # Search all subdirectories
  full.names = TRUE # Return full paths
)

# Create a function to read one file
read_alignment_stats = function(file) {
  # Extract barcode name from file path
  sample = str_extract(file, "barcode[0-9]+")
  
  # Read tsv
  read_tsv(
    file,
    comment = "#",
    col_names = FALSE, # No meaningful column names
    show_col_types = FALSE # Suppress noise
  ) %>%
    # Keep only rows where column 1 contains "SN" (stats rows)
    filter(X1 == "SN") %>%
    mutate(
      sample = sample,
      stat = str_remove(X2, ":$"),
      value = as.numeric(X3),
      .keep = "none" # Only keep these columns
    )
}

# Read all files and combine
alignment_stats = files %>%
  # map() applies a function once per element and returns a list
  map(read_alignment_stats) %>%
  # list_rbind() takes a list of dataframes/tibbles and stacks them row-wise
  list_rbind()

# Create percent mapped
alignment_stats_percent = alignment_stats %>%
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  ) %>%
  mutate(percent_mapped = (`reads mapped` / (`reads mapped` + `reads unmapped`)) * 100, 0) %>%
  select("sample", "percent_mapped")

# Transpose to be merged later
# Extract values to be used as column headers in the transposed dataframe
colnames = alignment_stats_percent[["sample"]]
# Transpose the dataset - makes rows (proteins) the columns and columns (samples) as rows
alignment_stats_percent = t(alignment_stats_percent[ ,-1])
# Assign previous values as column headers
colnames(alignment_stats_percent) = colnames
# rownames_to_column moves the row names to a new column
alignment_stats_percent = rownames_to_column((as.data.frame(alignment_stats_percent)), var="stat")

stat_order = c(
  "reads mapped",
  "reads unmapped",
  "percent_mapped",
  "average length",
  "maximum length",
  "average quality"
)

# Pivot wider and clean
alignment_stats_wide = alignment_stats %>%
  pivot_wider(
    names_from = "sample",
    values_from = "value"
  ) %>%
  # Bind by column name (not position)
  bind_rows(alignment_stats_percent) %>%
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
alignment_stats_formatted = alignment_stats_wide %>%
  mutate(
    across(-stat, ~ case_when(
      stat == "percent_mapped" ~ number(.x, accuracy=0.1),
      stat == "average quality" ~ number(.x, accuracy=0.1),
      TRUE ~ comma(.x)
    ))
  )

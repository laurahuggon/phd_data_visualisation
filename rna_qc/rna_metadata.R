# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl) # Read Excel files into R


# Experimental ------------------------------------------------------------

experimental_metadata = read_excel("/Users/k21224575/Library/CloudStorage/OneDrive-King\'sCollegeLondon/phd/lab/omics/rna/ipsc/rna_samples.xlsx", sheet = 2)

experimental_metadata = experimental_metadata %>%
  mutate(Fraction = case_when(
    Fraction == "Cytosol" ~ "Cyt",
    Fraction == "Synaptosome" ~ "Syn",
    TRUE ~ Fraction)) %>%
  mutate(Sample = paste(Genotype, Fraction, BiolRep, sep="_"))

experimental_metadata_collapsed = experimental_metadata %>%
  group_by(Sample, Genotype, Fraction, BiolRep, Working_stock, Date_plating, Date_differentiation, Date_RNA_extraction) %>%
  summarise(
    Avg_A260_A280 = round(mean(A260_A280), 2),
    Avg_A260_A230 = round(mean(A260_A230), 2),
    Avg_RINe = round(mean(RINe), 1),
    Avg_yield_ng = round(mean(Yield_ng), 0)
  ) %>%
  arrange(BiolRep, Fraction, desc(Genotype))
  

# Sequencing stats --------------------------------------------------------
# Load data ---------------------------------------------------------------

# Find all files
files = list.files(
  path = "/Volumes/prj/bcn_synaptopathy_in_ftd_als/synaptosome_long_read/data/nanostat",
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
  arrange(stat)

# Reformat from scientific notation
sequencing_stats_formatted = sequencing_stats_wide %>%
  mutate(
    across(-stat, ~ case_when(
      stat == "longest_read_(with_Q):1" ~ as.character(.x),
      stat == "median_qual" ~ number(as.numeric(.x), accuracy=0.1),
      TRUE ~ comma(as.numeric(.x))
    ))
  )

# Transpose
#' Transpose dataframe - old column headers are a column in the transposed dataframe
transpose_df_col = function(data, col_name, new_col_name, cols_to_exclude=c(1)){ # cols_to_exclude in the format c(x)
  # Extract values to be used as column headers in the transposed dataframe
  colnames = data[[col_name]]
  # Transpose the dataset - makes rows (proteins) the columns and columns (samples) as rows
  data = t(data[ ,-cols_to_exclude])
  # Assign previous values as column headers
  colnames(data) = colnames
  # rownames_to_column moves the row names to a new column
  data = rownames_to_column((as.data.frame(data)), var=new_col_name)
}
  
sequencing_stats_formatted = transpose_df_col(sequencing_stats_formatted, "stat", "Sample")

sequencing_stats_formatted = sequencing_stats_formatted %>%
  rename(
    Number_of_passed_reads = number_of_reads,
    Median_read_length = median_read_length,
    N50 = n50,
    Longest_read_with_Q = `longest_read_(with_Q):1`,
    Median_quality = median_qual
  )


# Alignment stats ---------------------------------------------------------
# Load data ---------------------------------------------------------------

# Find all files
files = list.files(
  path = "/Volumes/prj/bcn_synaptopathy_in_ftd_als/synaptosome_long_read/data/processed",
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
alignment_stats_percent = transpose_df_col(alignment_stats_percent, "sample", "stat")

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
  arrange(stat)

# Reformat from scientific notation
alignment_stats_formatted = alignment_stats_wide %>%
  mutate(
    across(-stat, ~ case_when(
      stat == "percent_mapped" ~ number(.x, accuracy=0.1),
      stat == "average quality" ~ number(.x, accuracy=0.1),
      TRUE ~ comma(.x)
    ))
  )

alignment_stats_formatted = transpose_df_col(alignment_stats_formatted, "stat", "Sample")

alignment_stats_formatted = alignment_stats_formatted %>%
  rename(
    Number_of_reads_mapped = `reads mapped`,
    Number_of_reads_unmapped = `reads unmapped`,
    Reads_aligned_to_genome = percent_mapped,
    Average_alignment_length = `average length`,
    Longest_alignment_length = `maximum length`,
    Average_quality = `average quality`
  )


# Merge -------------------------------------------------------------------

rna_sample_metadata = experimental_metadata_collapsed %>%
  left_join(sequencing_stats_formatted, by="Sample") %>%
  left_join(alignment_stats_formatted, by="Sample")

write_csv(rna_sample_metadata, "../synaptosome_long-read/metadata/rna_sample_metadata.csv")

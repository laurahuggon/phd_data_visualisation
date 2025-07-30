
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(scales)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King\'sCollegeLondon/phd/lab/omics/rna/ipsc/lrs/MinKNOW/barcodes/"

# Load data ---------------------------------------------------------------

# List all files in the directory
files = list.files(parent_filepath)
ids = str_extract(files, "PB[A-Za-z0-9]+")


# Read each file into a named list of dataframes
barcode_dfs = setNames(
  lapply(paste0(parent_filepath, files), read.csv, header = TRUE),
  ids
)

# Data processing ---------------------------------------------------------

# Function to calculate the number of reads passed
calculate_pass = function(df){
  df = df %>%
    mutate(total_passed_reads = (`Total.reads..k.`/100)*`Passed.reads....`)
  
  return(df)
}

# Apply to all dataframes
barcode_dfs = lapply(barcode_dfs, calculate_pass)
barcode_dfs[["PBC83182"]]

# Create one dataframe with number of reads passed for each sample
# Take the first column from the first dataframe
passed_reads = barcode_dfs[[1]][, 1, drop = FALSE]  # Keeps it as a data frame

# Extract total_passed_reads from each df and combine
reads_matrix = sapply(barcode_dfs, function(df) df$total_passed_reads)

# Bind them together
passed_reads = cbind(passed_reads, as.data.frame(reads_matrix))

# Only 16 samples
passed_reads = passed_reads[1:16, ]

# Sum across all flow cells and multiple by 1000
passed_reads = passed_reads %>%
  mutate(total = rowSums(across(-1))) %>% # across(-1) selects all columns except the first
  mutate(total = total*1000)

# Add sample columns
sample_names = c("WT_Cyt_1", "Q331K_Cyt_1", "WT_Syn_1", "Q331K_Syn_1", "WT_Cyt_2", "Q331K_Cyt_2", "WT_Syn_2", "Q331K_Syn_2", "WT_Cyt_3", "Q331K_Cyt_3", "WT_Syn_3", "Q331K_Syn_3", "WT_Cyt_4", "Q331K_Cyt_4", "WT_Syn_4", "Q331K_Syn_4")

passed_reads = passed_reads %>%
  mutate(sample = sample_names) %>%
  mutate(genotype = str_extract(sample, "^[^_]+")) %>%
  mutate(fraction = str_extract(sample, "(?<=_)[^_]+(?=_)"))

passed_reads$genotype = factor(passed_reads$genotype, levels=c("WT", "Q331K"))
passed_reads$sample = factor(passed_reads$sample, levels=c("WT_Cyt_1", "WT_Cyt_2", "WT_Cyt_3", "WT_Cyt_4", "Q331K_Cyt_1", "Q331K_Cyt_2", "Q331K_Cyt_3", "Q331K_Cyt_4", "WT_Syn_1", "WT_Syn_2", "WT_Syn_3", "WT_Syn_4", "Q331K_Syn_1", "Q331K_Syn_2", "Q331K_Syn_3", "Q331K_Syn_4"))

# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12), # Adjust plot title
          axis.title.x = element_text(margin = margin(t = 15), size = 12), # Adjust x-axis title
          axis.title.y = element_text(margin = margin(r = 15), size = 12), # Adjust y-axis title
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10), # Increase y-axis text size
          # Facet-specific
          panel.spacing = unit(0.5, "lines"), # Adjust spacing between facet panels
          strip.text = element_text(size = 10, face = "bold") # Facet title size
    ) 
}


# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, x, fill_group) {
  
  # Create the bar plot
  p = ggplot(group_data, aes(x = .data[[x]], 
                             y = .data[[group_col_name]],
                             fill = .data[[fill_group]])) +
    
    # Bar plot
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF", "Cyt" = "#A1D0E6", "Syn" = "#91ACA1")) +
    
    # Graph titles
    labs(title = "Sequencing depth per sample",
         x = "",
         y = "Passed reads",
         fill = "Genotype") + # Legend title

    # Plot appearance
    my_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 45000000), expand = c(0, 0), labels = label_number(scale = 1e-6, suffix = "M")) # Setting both multiplier and add-on to 0
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(passed_reads, "total", "sample", "fraction")
plot

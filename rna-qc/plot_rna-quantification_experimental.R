# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/omics/"
relative_filepath = "rna/ipsc/"
filename = "rna_samples.xlsx"
title = ""

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
nanodrop_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Find means --------------------------------------------------------

# # Calculate mean for each group
# nanodrop_means = nanodrop_data %>%
#   group_by(sample) %>%
#   summarise(
#     n = n(), # Count the number of samples in each group
#     concentration_mean = round(mean(concentration), 1),
#     a260_a280_mean = round(mean(`a260/a280`), 2),
#     a260_a230_mean = round(mean(`a260/a230`), 2),
#     RINe = RINe
#     ) 
# 
# # Split the sample column
# nanodrop_means = nanodrop_means %>%
#   separate(sample, into = c("fraction", "condition"), sep = " ", extra = "merge")

# Create sample_name column
nanodrop_data = nanodrop_data %>%
  mutate(Fraction = recode(Fraction,
                           "Cytosol" = "Cyt",
                           "Synaptosome" = "Syn")) %>%
  mutate(Sample_Name = paste0(Genotype, "_", Fraction, "_", BiolRep, ".", TechRep))

# Arrange
nanodrop_data = nanodrop_data %>%
  arrange(
    factor(Fraction, levels = c("Cyt", "Syn")),
    factor(Genotype, levels = c("WT", "Q331K"))
  ) %>%
  # Define as factor
  mutate(Sample_Name = factor(Sample_Name, levels = unique(Sample_Name)))

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
          strip.text = element_text(size = 10) # Facet title size
    ) 
}

# Create a function that takes one dataframe and column name to generate faceted bar plots by label
plot_data = function(group_data=nanodrop_data, y, x="Sample_Name", fill_group="Fraction") {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(nanodrop_data[[y]], na.rm = TRUE)
  upper_limit = max_y_value * 1.1 # 10% buffer above the max value
  
  # Compute y-axis title
  if(y == "ng_uL"){
    y_title = "Concentration (ng/μL)"
  } else if (y == "A260_A280"){
    y_title = "A260/A280"
  } else if (y == "A260_A230"){
    y_title = "A260/A230"
  } else if (y == "RINe"){
    y_title = "RINe"
  }
  
  # Create the bar plot
  p = ggplot(group_data, aes(x = .data[[x]], 
                             y = .data[[y]],
                             fill = .data[[fill_group]])) +
    
    # Bar plot
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF", "Cyt" = "#A1D0E6", "Syn" = "#91ACA1")) +
    
    # Graph titles
    labs(#title = title,
      x = "",
      y = y_title,
      fill = "Fraction") +
    
    # Plot appearance
    my_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add horizontal line at y = 1 and change y-axis
  if (y == "ng_uL"){
    p = p + scale_y_continuous(limits = c(0, upper_limit), breaks = seq(0, upper_limit, by = 50), expand = c(0, 0)
    )
  } else if (y == "A260_A280" || y == "A260_A230"){
    p = p + geom_hline(yintercept=1.8, linetype="dashed", color="black", size=0.3)
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))
  } else if (y == "RINe"){
    p = p + geom_hline(yintercept=8, linetype="dashed", color="black", size=0.3)
    p = p + scale_y_continuous(limits = c(0, upper_limit), breaks=seq(0, 10, 1), expand = c(0, 0))
  }
  
  # Print the plot
  return(p)
}

# Make plot
plot_conc = plot_data(y = "ng_uL")
plot_a280 = plot_data(y = "A260_A280")
plot_a230 = plot_data(y = "A260_A230")
plot_rin = plot_data(y = "RINe")

plot_conc
plot_a280
plot_a230
plot_rin

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, "rna_samples_plot_conc.png"), plot=plot_conc, width=18, height=4, dpi=300, bg="white")
# Save plot
ggsave(paste0(parent_filepath, relative_filepath, "rna_samples_plot_a280.png"), plot=plot_a280, width=18, height=4, dpi=300, bg="white")
# Save plot
ggsave(paste0(parent_filepath, relative_filepath, "rna_samples_plot_a230.png"), plot=plot_a230, width=18, height=4, dpi=300, bg="white")
# Save plot
ggsave(paste0(parent_filepath, relative_filepath, "rna_samples_plot_rin.png"), plot=plot_rin, width=18, height=4, dpi=300, bg="white")

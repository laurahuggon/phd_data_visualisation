
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/i3neuron_synapse/synper_extraction_practice/"
relative_filepath = "hmr/"
filename = "empiria_hmr.xlsx"
protein_name = "Homer-1"

colour1 = "#A1D0E6"
colour2 = "#91ACA1"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Prepare data frame
empiria_data = empiria_data[!is.na(empiria_data[, 4]), ] # Directly removes any rows with any `NA` values in the fourth column (Normalised Signal)
                                                         # Only keep the rows where fourth column is not NA
colnames(empiria_data) = as.character(empiria_data[1,]) # Make the first row the column names
empiria_data = empiria_data[-1, ] # Remove first row

# Split the Name column
empiria_data = empiria_data %>%
  separate(Name, into = c("Replicate", "Fraction"), sep = " ")

# Define `Normalized Signal` variable as numeric
empiria_data$`Normalized Signal` = as.numeric(empiria_data$`Normalized Signal`)

# Define Replicate and Fraction as a factor with levels
empiria_data$Replicate = factor(empiria_data$Replicate, levels = c("WT", "Q331K"))
empiria_data = empiria_data %>%
  mutate(Fraction = recode(Fraction,
                           "Homogenate" = "Hom",
                           "Cytosol" = "Cyt",
                           "Synaptosome" = "Syn"))
empiria_data$Fraction = factor(empiria_data$Fraction, levels = c("Hom", "Cyt", "Syn"))

# Remove homogenate
empiria_data = empiria_data %>%
  filter(Fraction != "Hom") %>%
  select(Lane, Replicate, Fraction, MW, `Normalized Signal`)


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for facet bar plots
my_theme_facet = function() {
  theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          axis.title.y = element_text(margin = margin(r = 15), # Adjust y-axis title position
                                      size = 12), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10), # Increase y-axis text size
          # Facet-specific
          panel.spacing = unit(1, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 11, face = "bold"),  # Increase facet title size and make it bold
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_normalised = function(data, x = "Fraction", y = "`Normalized Signal`") {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(empiria_data$`Normalized Signal`, na.rm = TRUE)
  upper_limit = max_y_value * 1.1 # 10% buffer above the max value
  
  # Create the bar plot
  p = ggplot(empiria_data, aes_string(x = x,
                                    y = y,
                                    fill = "Fraction")) +
    # Bar plot
    geom_col(position = position_dodge(0.9),
             width = 0.8,
             color = "black") +
    scale_fill_manual(values = c("Cyt" = colour1, "Syn" = colour2)) +
    
    # Facet by Replicate
    facet_wrap(~Replicate, axes="all") +
    
    # Graph titles
    labs(title = protein_name,
         x = "",
         y = "Normalised Signal",
         fill = x) +
    
    # Plot appearance
    my_theme_facet() +
    scale_y_continuous(limits = c(0, upper_limit),   # Set the y-axis range
                       expand = c(0, 0)  # Remove extra space around the axis
                       ) +
    
    # Overlay individual data points
    geom_point(data = empiria_data, aes_string(x = x,
                                               y = y),
                       position = position_dodge(0.9), size = 1.5)
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_normalised(empiria_data)

plot

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_enrichment_cytsyn.png"), plot=plot, width=3.2, height=3.5, dpi=300, bg="white")

# # Print message
# print(message)

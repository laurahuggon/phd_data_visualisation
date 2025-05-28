
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

# Rename column
empiria_data = empiria_data %>%
  rename(
    Normalised_Signal = `Normalized Signal`
  )


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12), # Adjust plot title
          axis.title.x = element_text(margin = margin(t = 15), size = 12), # Adjust x-axis title
          axis.title.y = element_text(margin = margin(r = 15), size = 12), # Adjust y-axis title
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10), # Increase y-axis text size
          # Facet-specific
          panel.spacing = unit(0.5, "lines"), # Adjust spacing between facet panels
          strip.text = element_text(size = 12, face = "bold") # Facet title size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, individual_data, individual_col_name, x, facet_grouping) {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(individual_data[[individual_col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.1  # 25% buffer above the max value
  
  # Reformulate
  facet_formula = reformulate(facet_grouping)
  
  # Create the bar plot
  p = ggplot(group_data, aes_string(x = x, 
                                    y = group_col_name,
                                    fill = x)) +
    
    # Bar plot
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("Cyt" = colour1, "Syn" = colour2)) +
    
    # Error bars
    # geom_errorbar(aes(ymin = .data[[group_col_name]] - SD,
    #                   ymax = .data[[group_col_name]] + SD),
    #               width = 0.2) +
    
    # Facet
    facet_wrap(facet_formula, nrow=1, axes= "all") +
    
    # Graph titles
    labs(title = protein_name,
         x = "",
         y = "Normalised Signal",
         fill = x) + # Legend title
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Define shapes for each DIFF value
  diff_shapes = c("6" = 21, "7" = 22, "8" = 24)
  
  # Overlay individual data points (optional with different shapes for DIFF)
  p = p + geom_point(data=individual_data, aes_string(x = x,
                                                            y = individual_col_name
                                                            ),
                           size = 1.25,
                           fill = "black")
    # scale_shape_manual(values = diff_shapes) 
  
  # # Significance stars
  # p = p + geom_text(data = test_results, aes(label = stars,
  #                                            x = 1.5,
  #                                            y = (max_y + 0.1*upper_limit)),
  #                   inherit.aes = FALSE,
  #                   size = 6)  # Adjust size here
  # # Significance lines
  # p = p + geom_segment(data = test_results, aes(x = 1,
  #                                               xend = 2,
  #                                               y = (max_y + 0.075*upper_limit),
  #                                               yend = (max_y + 0.075*upper_limit)),
  #                      linetype = "solid",
  #                      color = "black",
  #                      inherit.aes = FALSE) 	
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(empiria_data, "Normalised_Signal", empiria_data, "Normalised_Signal", "Fraction", "Replicate")

plot

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_enrichment_cytsyn.png"), plot=plot, width=3, height=3.5, dpi=300, bg="white")

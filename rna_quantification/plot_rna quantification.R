
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/rna/"
relative_filepath = "primaries/"
filename = "prep_2.xlsx"
title = "Optimisation 1 | Pooled 2x wells | Eluted in 10 μL"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
nanodrop_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Find means --------------------------------------------------------

# Calculate mean for each group
nanodrop_means = nanodrop_data %>%
  group_by(sample) %>%
  summarise(
    n = n(), # Count the number of samples in each group
    concentration_mean = round(mean(concentration), 1),
    a260_a280_mean = round(mean(`a260/a280`), 2),
    a260_a230_mean = round(mean(`a260/a230`), 2),
    RINe = RINe
    ) 

# Split the sample column
nanodrop_means = nanodrop_means %>%
  separate(sample, into = c("fraction", "condition"), sep = " ", extra = "merge")

# Define as factors
nanodrop_means$fraction = factor(nanodrop_means$fraction, levels = c("Homogenate", "Cytosol", "Synaptosome"))
nanodrop_means$condition = factor(nanodrop_means$condition, levels = c("-P -R", "+P -R", "+P +R", "-P +R"))


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for facet bar plots
my_theme_facet = function() {
  theme_minimal() +
    theme(
          plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          panel.spacing = unit(1, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 12, face = "bold"),  # Increase facet title size and make it bold
          axis.title.y = element_text(margin = margin(r = 10), # Adjust y-axis title position
                                      size = 12), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Create a function that takes one dataframe and column names to generate multiple bar plots with overlayed data points
plot = function(nanodrop_means, x = "fraction", y, fill = "condition") {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(nanodrop_means[[y]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25 # 25% buffer above the max value
  
  # Compute y-axis title
  if(y == "concentration_mean"){
    y_title = "Concentration (ng/μL)"
  } else if (y == "a260_a280_mean"){
    y_title = "A260/A280"
  } else if (y == "a260_a230_mean"){
    y_title = "A260/A230"
  } else if (y == "RINe"){
    y_title = "RINe"
  }
  
  # Create the bar plot
  p = ggplot(nanodrop_means, aes_string(x = x,
                                        y = y,
                                        fill = fill)) +
    # Bar plot
    geom_col(position = position_dodge(0.7),
             width = 0.6,
             color = "black") +
    
    # Graph titles
    labs(#title = title,
         x = "",
         y = y_title,
         fill = "Condition") +
    
    # Plot appearance
    my_theme_facet() +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Add horizontal line at y = 1
  if(y == "a260_280_mean" || y == "a260_a230_mean"){
    p = p + geom_hline(yintercept=1.8, linetype="dashed", color="black", size=0.3)
  } else if (y == "RINe"){
    p = p + geom_hline(yintercept=7, linetype="dashed", color="black", size=0.3)
    p = p + scale_y_continuous(limits = c(0, upper_limit), breaks=seq(2, 10, 2), expand = c(0, 0))
  }
  
  # Print the plot
  return(p)
}

# Make plot
plot_conc = plot(nanodrop_means, y = "concentration_mean")
plot_a280 = plot(nanodrop_means, y = "a260_a280_mean")
plot_a230 = plot(nanodrop_means, y = "a260_a230_mean")
plot_rin = plot(nanodrop_means, y = "RINe")

plot_conc
plot_a280
plot_a230
plot_rin


# Export plot
filename_noextn = str_split(filename, "\\.", simplify = TRUE)[1]

export = function(plot_name){
  # Create name  
  name = deparse(substitute(plot_name))
  
  # Open a PNG file to save the plot
  png(paste0(parent_filepath, relative_filepath, filename_noextn, "_", name, ".png"), width=1735, height=1335, res=300)
  
  # Create a plot
  print(plot_name)
  
  # Close the device
  dev.off()
}

export(plot_conc)
export(plot_a280)
export(plot_a230)
export(plot_rin)

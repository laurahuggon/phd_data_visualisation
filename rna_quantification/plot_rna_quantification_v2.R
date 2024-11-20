# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/omics/"
relative_filepath = "rna/ipsc/"
filename = "ipsc_characterisation.xlsx"
title = ""

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, filename)
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

# Define as factors
nanodrop_data$Genotype = factor(nanodrop_data$Genotype, levels = c("WT", "Q331K"))
nanodrop_data$Fraction = factor(nanodrop_data$Fraction, levels = c("Cytosol", "Synaptosome"), labels = c("Cyt", "Syn"))
nanodrop_data$Sample = factor(nanodrop_data$Sample)
# nanodrop_means$condition = factor(nanodrop_means$condition, levels = c("-P -R", "+P -R", "+P +R", "-P +R"))

# Create a new column for the label, combining 'Fraction' and 'BiolRep'
nanodrop_data$Label <- paste(nanodrop_data$Fraction, nanodrop_data$BiolRep, sep = " ")
nanodrop_data$Label = factor(nanodrop_data$Label, levels = c("Cyt 1", "Cyt 2", "Cyt 3", "Cyt 4", "Syn 1", "Syn 2", "Syn 3", "Syn 4"))

# Create a new column for the group, combining 'Genotype' and 'Fraction'
nanodrop_data$Group <- paste(nanodrop_data$Genotype, nanodrop_data$Fraction, sep = " ")
nanodrop_data$Group = factor(nanodrop_data$Group, levels = c("WT Cyt", "Q331K Cyt", "WT Syn", "Q331K Syn"))

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
      axis.text.y = element_text(size = 10), # Increase y-axis text size
    ) 
}

# Create a function that takes one dataframe and column names to generate multiple bar plots with overlayed data points
plot = function(nanodrop_data, x = "Sample", y, fill = "Group") {
  
  # Verify that the specified column exists in the data
  if (!y %in% colnames(nanodrop_data)) {
    stop(paste("Column", y, "not found in nanodrop_data"))
  }
  
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
  p = ggplot(nanodrop_data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    
    # Bar plot
    geom_col(position = position_dodge(0.7),
             width = 0.6,
             color = "black") +
    scale_fill_manual(values = c("WT Cyt" = "royalblue4", "Q331K Cyt" = "skyblue", "WT Syn" = "violetred3", "Q331K Syn" = "violet")) +
    
    # Graph titles
    labs(#title = title,
      x = "",
      y = y_title,
      fill = "Genotype") +
    
    # Plot appearance
    my_theme_facet() +
    facet_wrap(~Label, scales = "free_x", strip.position = "bottom", nrow = 1) +
    # Adjust spacing between facets
    # Adjust spacing between facets
    theme(
      panel.spacing = unit(0.25, "lines"),  # Reduce the space between the panels
      strip.text.x = element_text(face = "plain"),  # Edit facet titles
      strip.placement = "outside",  # Moves the facet labels outside the plot area
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.ticks.x = element_blank()  # Remove x-axis ticks
    )
  
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
plot_conc = plot(nanodrop_data, y = "ng_uL")
plot_a280 = plot(nanodrop_data, y = "A260_A280")
plot_a230 = plot(nanodrop_data, y = "A260_A230")
plot_rin = plot(nanodrop_data, y = "RINe")

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
  png(paste0(parent_filepath, relative_filepath, filename_noextn, "_", name, ".png"), width=3500, height=1335, res=300)
  
  # Create a plot
  print(plot_name)
  
  # Close the device
  dev.off()
}

export(plot_conc)
export(plot_a280)
export(plot_a230)
export(plot_rin)
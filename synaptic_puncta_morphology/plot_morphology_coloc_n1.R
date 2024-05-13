
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y2/unc13/analysis_nis_elements/"
relative_filepath = "unc13_coloc/"
filename = "COLOCPRE_BTUB_1-2.csv"
entity = "PRE"
pre_marker = "UNC13A"
post_marker = "PSD-95"


# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
nis_elements_df = read_csv(full_filename)


# Prepare data ------------------------------------------------------------

# Unblind samples by adding Genotype, DIFF and DIV columns
# Define mappings - map each suffix to its corresponding value
if (pre_marker == "UNC13A" && post_marker == "PSD-95") {
  genotype_map = c("1" = "WT", "2" = "Q331K")
  
  diff_map = c("1" = "14", "2" = "14")
  
  div_map = c("1" = "21", "2" = "21")
  
  suffix_regex = "(?<=_)[1-2](?=_)" # Uses lookbehind `(?<=_)` and lookahead `(?=_)` to capture character between two underscores
  
}

# Extract the suffix from the filename
nis_elements_df = nis_elements_df %>%
  mutate(Suffix = str_extract(Filename, suffix_regex)) %>%
  mutate(
    Genotype = genotype_map[Suffix], # Map suffixes to respective values using predefined mappings
    DIFF = diff_map[Suffix],
    DIV = div_map[Suffix]
  ) %>%
  select(-Suffix)  # Remove the Suffix column if it's not needed later

# Remove DIV 28
nis_elements_df = nis_elements_df %>%
  filter(DIV != "28")

# Define Genotype and DIV variable as a factor with levels
nis_elements_df$Genotype = factor(nis_elements_df$Genotype, levels = c("WT", "Q331K"))
nis_elements_df$DIV = factor(nis_elements_df$DIV, levels = c("7", "14", "21"))


# Find sample means -------------------------------------------------------

# Create function that finds sample means for a given variable
mean_by_sample = function(data, column_name) {
  # Group by Genotype, DIV, DIFF
  result = data %>%
    group_by(Genotype, DIV, DIFF) %>%
    summarise(
      N = n(), # Count the number of images in each sample
      N_Mean = mean(.data[[column_name]]), # Calculate mean for each sample
    )
  return(result)
}

# Find mean colocalisation, density and volume for each sample
sample_mean_coloc = mean_by_sample(nis_elements_df, "Coloc")
sample_mean_density = mean_by_sample(nis_elements_df, "DensityColoc")
sample_mean_volume = mean_by_sample(nis_elements_df, "MeanVolumeColoc")

# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, column_name = "N_Mean") {
  # Group by Genotype, DIV
  result = data %>%
    group_by(Genotype, DIV) %>%
    summarise(
      N = n(), # Count the number of samples in each group
      Global_Mean = mean(.data[[column_name]]), # Calculate mean for each group
      SD = sd(.data[[column_name]]) # Calculate the SD for each group
    )
  return(result)
}

# Find mean colocalisation, density and volume for each group
group_mean_coloc = mean_by_group(sample_mean_coloc)
group_mean_density = mean_by_group(sample_mean_density)
group_mean_volume = mean_by_group(sample_mean_volume)


# Create a function that calculates the maximum y-values per facet -> this is for dynamic annotation bars in the plots
calculate_max_y_per_facet = function(group_data) {
  # Add a new column to the group_data that calculates the potential max height for the error bars
  group_data$max_y = group_data$Global_Mean
  
  # Aggregate these max heights by DIV to get the maximum for each facet
  max_y_per_div = aggregate(max_y ~ DIV, data = group_data, max)
  
  return(max_y_per_div)
}

# Find maximum y-values per facet for each data type
max_y_per_div_coloc = calculate_max_y_per_facet(group_mean_coloc)
max_y_per_div_density = calculate_max_y_per_facet(group_mean_density)
max_y_per_div_volume = calculate_max_y_per_facet(group_mean_volume)


# # Statistics --------------------------------------------------------------
# 
# # Create function to perform unpaired two-sided t-tests for each DIV value
# perform_t_tests = function(data) {
#   # Extract unique DIV levels to perform the t-test for each level separately
#   div_levels = unique(data$DIV) # Extract unique values from the DIV column
#   
#   # Create an empty list to store test results
#   test_results = list()
#   
#   # Loop through each DIV level - allows performance of separate analyses at each time point
#   for (div in div_levels) {
#     # Subset data by DIV
#     data_div = subset(data, DIV == div) # Extract rows from dataframe where DIV matches the current DIV level in the loop
#     
#     # Extract N_Mean values for each genotype
#     wt_data = subset(data_div, Genotype == "WT")$N_Mean
#     q331k_data = subset(data_div, Genotype == "Q331K")$N_Mean
#     
#     # Perform Student's t-test
#     # Function call:
#     # `wt_data` and `q331k_data` are vectors containing the `N_Mean` values for each genotype at the current DIV level
#     # 
#     # Parameters:
#     # `var.equal = TRUE` assumes equal variances of the two groups
#     # `paired = FALSE` is the default
#     #
#     # Assignment to `test_results`:
#     # Store the results in a list, `test_results`, under a key created using the current DIV value (`div`)
#     test_results[[div]] = t.test(wt_data, q331k_data, alternative = "two.sided", var.equal = TRUE)
#   }
#   
#   # Return the test results
#   return(test_results)
# }
# 
# # Perform t-test for each data type
# coloc_test_results = perform_t_tests(sample_mean_coloc)
# density_test_results = perform_t_tests(sample_mean_density)
# volume_test_results = perform_t_tests(sample_mean_volume)
# 
# # Create function that converts p-values to signficance stars and generates a dataframe of annotations
# prepare_annotations = function(test_results) {
#   # Convert p-values to stars based on traditional significance levels
#   convert_p_to_stars = function(p_value) {
#     if (p_value <= 0.0001) {
#       return("****")
#     } else if (p_value <= 0.001) {
#       return("***")
#     } else if (p_value <= 0.01) {
#       return("**")
#     } else if (p_value <= 0.05) {
#       return("*")
#     } else {
#       return("")  # Not significant
#     }
#   }
#   
#   # Create a dataframe `annotations`
#   annotations = data.frame(
#     DIV = names(test_results), # One entry per DIV level
#     p_value = sapply(test_results, function(x) x$p.value), # Extract p-values
#     stringsAsFactors = FALSE
#   )
#   
#   # Convert p-values to stars
#   annotations$Stars = sapply(annotations$p_value, convert_p_to_stars)
#   
#   # Set factor levels for `DIV`
#   annotations$DIV = factor(annotations$DIV, levels = c("7", "14", "21"))
#   
#   # Filter out non-significant annotations
#   annotations = annotations[annotations$Stars != "", ]
#   
#   return(annotations)
# }
# 
# # Create annotations for each data type
# coloc_annotations = prepare_annotations(coloc_test_results)
# density_annotations = prepare_annotations(density_test_results)
# volume_annotations = prepare_annotations(volume_test_results)
# 
# # Create a function that merges max y-values with annotations -> this is for dynamic annotation bars in the plots
# merge_annotations_with_max_y = function(annotations, max_y_per_div) {
#   # Merge the max y-values per DIV with the annotations data frame
#   annotations = merge(annotations, max_y_per_div, by = "DIV")
#   
#   return(annotations)
# }
# 
# # Merge max y-values with annotations for each data type
# coloc_annotations = merge_annotations_with_max_y(coloc_annotations, max_y_per_div_coloc)
# density_annotations = merge_annotations_with_max_y(density_annotations, max_y_per_div_density)
# volume_annotations = merge_annotations_with_max_y(volume_annotations, max_y_per_div_volume)


# Data visualisation ------------------------------------------------------

# Define the protein name for y-axis label
protein_name = if (entity == "PRE") {
  pre_marker
} else {
  post_marker
}

# Create custom ggplot2 theme for facet bar plots
my_theme_facet = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          panel.spacing = unit(1, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 12, face = "bold"),  # Increase facet title size and make it bold
          axis.title.y = element_text(margin = margin(r = 13.5), # Adjust y-axis title position
                                      size = 12), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_by_genotype_div = function(group_data, sample_data, x = "Genotype", sd = "SD", facet = "DIV") {
  # Extracting the name of the dataframe
  data_name = deparse(substitute(group_data)) # Get the name of the `group_data` dataframe as a string
  
  # Determine y-axis label based on the inferred dataframe name
  if (grepl("coloc", data_name)) {
    y_label = bquote("Colocalised " * .(protein_name) * " puncta (%)")
  } else if (grepl("density", data_name)) {
    y_label = bquote("Density of coloc " * .(protein_name) * " puncta (puncta/μm3)")
  } else if (grepl("volume", data_name)) {
    y_label = bquote("Volume of coloc " * .(protein_name) * " puncta (μm3)")
  } else {
    y_label = "Mean"  # Default label if no specific identifier is found
  }
  
  # Check if the necessary columns exist in the group_data
  if (!all(c(x, y = "Global_Mean", sd, facet) %in% names(group_data))) {
    stop("group_data does not contain the necessary columns.")
  }
  
  # Check if the necessary columns exist in the sample_data
  if (!all(c(x, y = "N_Mean", facet) %in% names(sample_data))) {
    stop("sample_data does not contain the necessary columns.")
  }
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(group_data[["Global_Mean"]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Define custom labeller function to add "DIV " to facet titles
  my_labeller = as_labeller(function(x) paste("DIV", x))
  
  # Create the bar plot
  p = ggplot(group_data, aes_string(x = x,
                                    y = "Global_Mean",
                                    fill = x)) +
    # Bar plot
    geom_col(position = position_dodge(0.9),
             width = 0.6,
             color = "black") +
    scale_fill_manual(values = c("WT" = "grey40", "Q331K" = "grey88")) +
    
    # Facet
    facet_wrap(as.formula(paste0("~ ", facet)), ncol = 3, labeller = my_labeller) +
    
    # Graph titles
    labs(x = "",
         y = y_label,
         fill = x) +
    
    # Plot appearance
    my_theme_facet() +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Overlay individual data points
  p = p + geom_point(data = sample_data, aes_string(x = x,
                                                    y = "N_Mean"),
                     position = position_dodge(0.9), size = 1.5)
  
  # Print the plot
  return(p)
}

# Make plots for each data type
coloc_plot = plot_by_genotype_div(group_mean_coloc, sample_mean_coloc)
density_plot = plot_by_genotype_div(group_mean_density, sample_mean_density)
volume_plot = plot_by_genotype_div(group_mean_volume, sample_mean_volume)

# Arrange plots vertically
# all_plots = coloc_plot /
#             density_plot /
#             volume_plot

# Arrange plots horizontally
all_plots = coloc_plot + density_plot + volume_plot

# Export plots
# Open a PNG file to save the plot

# For vertical arrangment
# png(paste0(filepath, protein_name, "_all_plots.png"), width=820, height=3800, res=300)

# For horizontal arrangement
png(paste0(parent_filepath, relative_filepath, protein_name, "_coloc_all_plots.png"), width=2452, height=1335, res=300)

# Create a plot
all_plots

# Close the device
dev.off()
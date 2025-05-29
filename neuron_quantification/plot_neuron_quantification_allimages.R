
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggbeeswarm)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King\'sCollegeLondon/phd/lab/imaging/echo/imaging_data_y1/neuron_quantification/neun/qupath/qupath_project/"


# Load data ---------------------------------------------------------------

# DAPI (Live)+ cells
dapi_files = list.files(paste0(parent_filepath, "annotation_tables_blue"))
dapi_df = data.frame()
for (file in dapi_files) {
  # Create filepath
  filepath = paste0(parent_filepath, "annotation_tables_blue/", file)
  # Read file and append to the dataframe
  dapi_df = rbind(dapi_df, read.table(filepath, header=TRUE, sep="\t"))
}

# DAPI (Live)+ BTUB+ cells
dapibtub_files = list.files(paste0(parent_filepath, "annotation_tables_redblue"))
dapibtub_df = data.frame()
for (file in dapibtub_files) {
  # Create filepath
  filepath = paste0(parent_filepath, "annotation_tables_redblue/", file)
  # Read file and append to the dataframe
  dapibtub_df = rbind(dapibtub_df, read.table(filepath, header=TRUE, sep="\t"))
}

# Merge dataframes
dapibtub_df = merge(dapi_df, dapibtub_df, by="Image")


# Prepare data ------------------------------------------------------------

# Extract DIFF
dapibtub_df$DIFF = as.numeric(str_extract(dapibtub_df$Image, "(?<=_image_)\\d+"))

# Extract genotype
dapibtub_df$Genotype = str_extract(dapibtub_df$Image, "(?<=_image_\\d_)(QK|WT)")
# Replace QK with Q331K
dapibtub_df$Genotype = gsub("QK", "Q331K", dapibtub_df$Genotype)

# Define Genotype and DIFF variable as a factor with levels
dapibtub_df$Genotype = factor(dapibtub_df$Genotype, levels = c("WT", "Q331K"))
dapibtub_df$DIFF = factor(dapibtub_df$DIFF, levels = c(6, 7, 8))

# Find percentage + cells
dapibtub_df = dapibtub_df %>%
  mutate(
    percentage_positive = (Num.Red..Blue / Num.Blue) * 100
  )


# Find sample means -------------------------------------------------------

# # Create function that finds sample means for a given variable
# mean_by_sample = function(data, value_col, grouping_col) {
#   # Group and summarise
#   result = data %>%
#     group_by(!!!syms(grouping_col)) %>%
#     summarise(
#       N = n(), # Count the number of images in each sample
#       Sample_Mean = mean(.data[[value_col]]), # Calculate mean for each sample
#     )
#   return(result)
# }
# 
# # Find mean for each sample
# sample_means = mean_by_sample(dapibtub_df, "percentage_positive", c("Genotype", "DIFF"))


# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, value_col, grouping_col) {
  # Group and summarise
  result = data %>%
    group_by(!!!syms(grouping_col)) %>%
    summarise(
      N = n(), # Count the number of samples in each group
      Group_Mean = mean(.data[[value_col]]), # Calculate mean for each group
      SD = sd(.data[[value_col]]) # Calculate the SD for each group
    )
  return(result)
}

# Find mean for each group
group_means = mean_by_group(dapibtub_df, "percentage_positive", "Genotype")


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_test = function(data, value_col, grouping_col, ctrl_group, exp_group) {
  
  # Extract data for groups
  # [[1]] pulls the first (and only) column as a vector
  ctrl_data = data[data[[grouping_col]] == ctrl_group, value_col]
  exp_data = data[data[[grouping_col]] == exp_group, value_col]
  
  # Test for normality using Shapiro-Wilk test
  # The variable is set to TRUE if the p-value is greater than 0.05 (if data is greater than 0.05, it is considered normally distributed)
  ctrl_normal = shapiro.test(ctrl_data)$p.value > 0.05
  exp_normal = shapiro.test(exp_data)$p.value > 0.05
  
  if (!ctrl_normal || !exp_normal) { # If either variable is set to `FALSE`
    # If data is not normally distributed, perform Mann-Whitney test
    test_result = wilcox.test(exp_data, ctrl_data)
  } else {
    # If both groups are normally distributed, test for equal variances using F-test
    # If p-value is greater than 0.05, it assumes both groups have equal variances
    var_test = var.test(exp_data, ctrl_data)$p.value > 0.05
    
    if (!var_test) { # If `var_equal` variable is set to `FALSE`
      # If variances are not equal, perform t-test with var.equal = FALSE
      test_result = t.test(exp_data, ctrl_data, var.equal = FALSE)
    } else {
      # If variances are equal, perform t-test with var.equal = TRUE					
      test_result = t.test(exp_data, ctrl_data, var.equal = TRUE)
    }
  }
  
  # Return test results
  return(list(test_result = test_result,
              n_ctrl = length(ctrl_data),
              n_exp = length(exp_data)))
}

# Perform t-test
test_result_list = perform_test(dapibtub_df, "percentage_positive", "Genotype", "WT", "Q331K")

# Convert to dataframe
test_result = data.frame(
  n_ctrl = test_result_list$n_ctrl,
  n_exp = test_result_list$n_exp,
  p.value = test_result_list$test_result$p.value,
  method = test_result_list$test_result$method,
  alternative = test_result_list$test_result$alternative
)

# Add significance stars
test_result = test_result %>%
  mutate(
    stars = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Find max y-value from individual data points
test_result = test_result %>%
  mutate(max_y = max(dapibtub_df$percentage_positive)) %>%
  # Replace max_y with NA if not significant
  mutate(max_y = ifelse(stars == "", NA, max_y))

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
          strip.text = element_text(12, face = "bold") # Facet title size
    ) 
} 

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, individual_data, individual_col_name, x, test_results = test_result) {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(individual_data[[individual_col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Create the bar plot
  p = ggplot(group_data, aes_string(x = x, 
                                    y = group_col_name,
                                    fill = x)) +
    
    # Bar plot
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF")) +
    
    # Error bars
    geom_errorbar(aes(ymin = .data[[group_col_name]] - SD,
                      ymax = .data[[group_col_name]] + SD),
                  width = 0.2) +
    
    # Graph titles
    labs(title = "",
         x = "",
         y = "β-III tubulin+ cells (%)",
         fill = x) + # Legend title
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Define shapes for each DIFF value
  diff_shapes = c("6" = 21, "7" = 22, "8" = 24)
  
  # Overlay individual data points (optional with different shapes for DIFF)
  p = p + geom_quasirandom(data=individual_data, aes_string(x = x,
                                                      y = individual_col_name,
                                                      shape = "DIFF"),
                           width = 0.2,
                     size = 1.25,
                     fill = "black", alpha = 0.5) +
    scale_shape_manual(values = diff_shapes) 
  
  # Significance stars
  p = p + geom_text(data = test_results, aes(label = stars,
                                             x = 1.5,
                                             y = (max_y + 0.1*upper_limit)),
                    inherit.aes = FALSE,
                    size = 6)  # Adjust size here
  # Significance lines
  p = p + geom_segment(data = test_results, aes(x = 1,
                                                xend = 2,
                                                y = (max_y + 0.075*upper_limit),
                                                yend = (max_y + 0.075*upper_limit)),
                       linetype = "solid",
                       color = "black",
                       inherit.aes = FALSE)
  
  # Print the plot
  return(p)
}


# Make plot
plot = plot_data(group_means, "Group_Mean", dapibtub_df, "percentage_positive", "Genotype")

plot

# Save plot
ggsave(paste0(parent_filepath, "btub_percentage_allimages.png"), plot=plot, width=1.85, height=3.5, dpi=300, bg="white")

# Export test results
# Define file path for saving CSVs
csv_path = paste0(parent_filepath, "test_result_allimages.csv")

# Export the test results to CSV files
write.csv(test_result, csv_path, row.names=FALSE)



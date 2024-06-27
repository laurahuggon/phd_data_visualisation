
# For normalising intensities to BTUB intensity: uncomment the calculation code, add "Normalised" the y-axis title on the plot, and change the file name from "_raw" to "_normalised".

# Load libraries ----------------------------------------------------------

library(tidyverse)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y2/unc13/analysis_nis_elements/cellbody_intensity/"

marker = "UNC13A"
entity = "PRE"
measurement = "SumIntensity"


# Load data ---------------------------------------------------------------

# Create relative_filepath and filename
if (marker == "UNC13A") {
  relative_filepath = "unc/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 3), ".csv")
} else if (marker == "synaptophysin" | marker == "PSD-95") {
  relative_filepath = "highdensity/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 11), ".csv")
}

full_filename = paste0(parent_filepath, relative_filepath, filename)
nis_elements_df = read_csv(full_filename)


# Prepare data ------------------------------------------------------------

# Unblind samples by adding Genotype, DIFF and DIV columns
# Define mappings - map each suffix to its corresponding value
if (marker == "UNC13A") {
  genotype_map = c("1" = "WT", "2" = "Q331K")
  
  diff_map = c("1" = "14", "2" = "14")
  
  div_map = c("1" = "21", "2" = "21")
  
  suffix_regex = "(?<=_)[1-2](?=_)" # Uses lookbehind `(?<=_)` and lookahead `(?=_)` to capture character between two underscores
} else if (marker == "synaptophysin" | marker == "PSD-95") {
  genotype_map = c("WT" = "WT", "QK" = "Q331K")
  
  diff_map = c("WT" = "3", "QK" = "3")
  
  div_map = c("WT" = "21", "QK" = "21")
  
  suffix_regex = "(?<=_)(WT|QK)(?=_)" # Uses lookbehind `(?<=_)` and lookahead `(?=_)` to capture either "WT" or "QK" between two underscores
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

# Define Genotype and DIFF variable as a factor with levels
nis_elements_df$Genotype = factor(nis_elements_df$Genotype, levels = c("WT", "Q331K"))
nis_elements_df$DIFF = factor(nis_elements_df$DIFF, levels = c(3, 4, 5, 14))

# Remove rows that found no cell body
nis_elements_df = nis_elements_df %>%
  filter(CellBodyCount != 0)

# Remove DIFF 5
if (marker == "syntaxin 1A" | marker == "Homer") {
  nis_elements_df = nis_elements_df %>%
    filter(DIFF != 5)
}

# Normalise SumIntensity values to TotalCellBodyVolume
if (entity == "PRE") {
  nis_elements_df = nis_elements_df %>%
    mutate(
      PRE_SumIntensity = PRE_SumIntensity / `TotalCellBodyVolume(BTUB)`
    )
} else if (entity == "POST") {
  nis_elements_df = nis_elements_df %>%
    mutate(
      POST_SumIntensity = POST_SumIntensity / `TotalCellBodyVolume(BTUB)`
    )
}

# Normalise MedianIntensity to BTUB_MedianIntensity
# if (entity == "PRE") {
#   nis_elements_df = nis_elements_df %>%
#     mutate(
#       PRE_MedianIntensity = PRE_MedianIntensity / BTUB_MedianIntensity
#     )
# } else if (entity == "POST") {
#   nis_elements_df = nis_elements_df %>%
#     mutate(
#       POST_MedianIntensity = POST_MedianIntensity / BTUB_MedianIntensity
#     )
# }

# Find sample means -------------------------------------------------------

# Create function that finds sample means for a given variable
mean_by_sample = function(data, column_name) {
  # Group by Genotype, DIFF
  result = data %>%
    group_by(Genotype, DIFF) %>%
    summarise(
      N = n(), # Count the number of images in each sample
      N_Mean = mean(.data[[column_name]]), # Calculate mean for each sample
    )
  return(result)
}

# Define column name
column_name = paste0(entity, "_", measurement)

# Find mean for each sample
sample_means = mean_by_sample(nis_elements_df, column_name)


# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, column_name = "N_Mean") {
  # Group by Genotype
  result = data %>%
    group_by(Genotype) %>%
    summarise(
      N = n(), # Count the number of samples in each group
      Global_Mean = mean(.data[[column_name]]), # Calculate mean for each group
      SD = sd(.data[[column_name]]) # Calculate the SD for each group
    )
  return(result)
}

# Find mean for each group
group_means = mean_by_group(sample_means)


# # Statistics --------------------------------------------------------------
# 
# # Create function to test normality, equal variance, and perform appropriate statistical test
# perform_test = function(data) {
#   
#   # Extract values grouped by Genotype using N_Mean
#   # Test for normality using Shapiro-Wilk test
#   normality_result = by(data$N_Mean, data$Genotype, shapiro.test)
#   
#   wt_normal = normality_result$WT$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
#   # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
#   q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
#   # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
#   
#   if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
#     # If data is not normally distributed, perform Mann-Whitney test
#     message = "Data not normally distributed; perform Mann-Whitney test"
#     test_result = wilcox.test(N_Mean ~ Genotype, data = data) # `N_Mean ~ Genotype` specifies that you want to compare values in the `N_Mean` column grouped by `Genotype`
#   } else {
#     # Test for equal variances using F-test
#     variance_test = var.test(N_Mean ~ Genotype, data = data)
#     var_equal = variance_test$p.value > 0.05 # If p-value is greater than 0.05, it assumes both genotypes have equal variances
#     
#     if (!var_equal) { # If `var_equal` variable is set to `FALSE`
#       # If variances are not equal, perform t-test with var.equal = FALSE
#       message = "Data does not have equal variance; perform Welch's t-test"
#       test_result = t.test(N_Mean ~ Genotype, data = data, var.equal = FALSE)
#     } else {
#       # If variances are equal, perform t-test with var.equal = TRUE
#       message = "Data has equal variance; perform Student's t-test"
#       test_result = t.test(N_Mean ~ Genotype, data = data, var.equal = TRUE)
#     }
#   }
#   
#   # Return both test results and the message
#   return(list(test_result = test_result, message = message))
# }
# 
# # Perform t-test
# result = perform_test(sample_means)
# 
# # Accessing the test results
# test_result = result$test_result
# 
# # Accessing the message
# message = result$message
# 
# # Create function that converts p-values to signficance stars and generates a dataframe of annotations
# prepare_annotations = function(test_result) {
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
#   # Extract p-value from test results
#   p_value = test_result$p.value
#   
#   # Convert p-value to stars
#   stars = convert_p_to_stars(p_value)
#   
#   # Create a dataframe `annotations`
#   annotations = data.frame(
#     p_value = p_value,
#     Stars = stars,
#     stringsAsFactors = FALSE
#   )
#   
#   # Filter out non-significant annotations
#   annotations = annotations[annotations$Stars != "", ]
#   
#   return(annotations)
# }
# 
# # Create annotations
# annotation = prepare_annotations(test_result)
# 
# # Find the maximum y-values -> this is for dynamic annotation bars in the plots
# annotation = annotation %>%
#   mutate(
#     max_y = max(group_means$Global_Mean + group_means$SD)
#   )


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold",
                                    hjust = 0.5), # Adjust plot title
          axis.title.y = element_text(margin = margin(r = 7.5), # Adjust y-axis title position
                                      size = 13), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, sample_data, x = "Genotype", sd = "SD") {
  # Check if the necessary columns exist in the group_data
  if (!all(c(x, y = "Global_Mean", sd) %in% names(group_data))) {
    stop("group_data does not contain the necessary columns.")
  }
  
  # Check if the necessary columns exist in the sample_data
  if (!all(c(x, y = "N_Mean") %in% names(sample_data))) {
    stop("sample_data does not contain the necessary columns.")
  }
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(group_data["Global_Mean"], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Create the bar plot
  p = ggplot(group_data, aes_string(x = x,
                                    y = "Global_Mean",
                                    fill = x)) +
    # Bar plot
    geom_col(position = position_dodge(0.9),
             width = 0.6,
             color = "black") +
    scale_fill_manual(values = c("WT" = "grey40", "Q331K" = "grey88")) +
    
    # Graph titles
    labs(title = "Cell Body",
         x = "",
         y = paste0(measurement, " of ", marker, " (a.u.)"),
         fill = x) +
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Define shapes for each DIFF value
  diff_shapes <- c("3" = 21, "4" = 22, "5" = 24, "14" = 25)
  
  # Overlay individual data points with different shapes for DIFF
  p = p + geom_point(data = sample_data, aes_string(x = x,
                                                    y = "N_Mean",
                                                    shape = "DIFF"),  # Added shape aesthetic
                     position = position_dodge(0), size = 1.5, fill = "black") +
    scale_shape_manual(values = diff_shapes)  # Adjust shape values if needed
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(group_means, sample_means)

plot

# Export plot
# Open a PNG file to save the plot
#
# For plot title with 1 line:
# width=825, height=1335
#
# For plot title with 2 lines:
# width=825, height=1390
png(paste0(parent_filepath, relative_filepath, marker, "_cellbody_", measurement, "_raw.png"), width=825, height=1335, res=300)

# Create a plot
plot

# Close the device
dev.off()

# Histograms
# plot2 = ggplot(nis_elements_df, aes(x = PRE_MeanIntensity, fill = DIFF, color = DIFF)) +
#   geom_density(alpha = 0.4, size = 0.5, adjust = 1.5) +  # Adjust the smoothness
#   labs(title = paste0("Density Plot of ", measurement, " (", marker, ")"),
#        x = measurement,
#        y = "Density") +
#   theme_minimal() +
#   theme(legend.position = "right") +
#   xlim(min(nis_elements_df$PRE_MeanIntensity) - 500, max(nis_elements_df$PRE_MeanIntensity) + 750) +
#   facet_wrap(~ Genotype, scales = "free_y")  # Separate plots by Genotype
# 
# plot2

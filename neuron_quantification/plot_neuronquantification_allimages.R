
# For normalising intensities to BTUB intensity: uncomment the calculation code, add "Normalised" the y-axis title on the plot, and change the file name from "_raw" to "_normalised".

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggbeeswarm)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King\'sCollegeLondon/phd/lab/imaging/echo/imaging_data_y1/neuron_quantification/neun/qupath/qupath_project/"


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
dapibtub_df$Diff = as.numeric(str_extract(dapibtub_df$Image, "(?<=_image_)\\d+"))

# Extract genotype
dapibtub_df$Genotype = str_extract(dapibtub_df$Image, "(?<=_image_\\d_)(QK|WT)")
# Replace QK with Q331K
dapibtub_df$Genotype = gsub("QK", "Q331K", dapibtub_df$Genotype)

# Define Genotype and DIFF variable as a factor with levels
dapibtub_df$Genotype = factor(dapibtub_df$Genotype, levels = c("WT", "Q331K"))
dapibtub_df$DIFF = factor(dapibtub_df$Diff, levels = c(6, 7, 8))

# Find percentage + cells
dapibtub_df = dapibtub_df %>%
  mutate(
    percentage_positive = (Num.Red..Blue / Num.Blue) * 100
  )

# Find sample means -------------------------------------------------------

# Create function that finds sample means for a given variable
# mean_by_sample = function(data, column_name) {
#   # Group by Genotype, DIFF
#   result = data %>%
#     group_by(Genotype, DIFF) %>%
#     summarise(
#       N = n(), # Count the number of images in each sample
#       N_Mean = mean(.data[[column_name]]), # Calculate mean for each sample
#     )
#   return(result)
# }
# 
# # Define column name
# column_name = paste0(entity, "_", measurement)
# 
# # Find mean for each sample
# sample_means = mean_by_sample(dapibtub_df, column_name)


# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, column_name) {
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
group_means = mean_by_group(dapibtub_df, "percentage_positive")


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_test = function(data, col_name) {
  
  # Extract values grouped by Genotype using col_name
  # Test for normality using Shapiro-Wilk test
  normality_result = by(data[[col_name]], data$Genotype, shapiro.test)
  
  wt_normal = normality_result$WT$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
  # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
  q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
  # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
  
  if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
    # If data is not normally distributed, perform Mann-Whitney test
    message = "Data not normally distributed; perform Mann-Whitney test"
    test_result = wilcox.test(data[[col_name]] ~ data$Genotype) # `col_name ~ Genotype` specifies that you want to compare values in the `col_name` column grouped by `Genotype`
  } else {
    # Test for equal variances using F-test
    variance_test = var.test(data[[col_name]] ~ data$Genotype)
    var_equal = variance_test$p.value > 0.05 # If p-value is greater than 0.05, it assumes both genotypes have equal variances
    
    if (!var_equal) { # If `var_equal` variable is set to `FALSE`
      # If variances are not equal, perform t-test with var.equal = FALSE
      message = "Data does not have equal variance; perform Welch's t-test"
      test_result = t.test(data[[col_name]] ~ data$Genotype, var.equal = FALSE)
    } else {
      # If variances are equal, perform t-test with var.equal = TRUE
      message = "Data has equal variance; perform Student's t-test"
      test_result = t.test(data[[col_name]] ~ data$Genotype, var.equal = TRUE)
    }
  }
  
  # Return both test results and the message
  return(list(test_result = test_result, message = message))
}

# Perform t-test
result = perform_test(dapibtub_df, "percentage_positive")

# Accessing the test results
test_result = result$test_result

# Accessing the message
message = result$message

# Create function that converts p-values to signficance stars and generates a dataframe of annotations
prepare_annotations = function(test_result) {
  # Convert p-values to stars based on traditional significance levels
  convert_p_to_stars = function(p_value) {
    if (p_value <= 0.0001) {
      return("****")
    } else if (p_value <= 0.001) {
      return("***")
    } else if (p_value <= 0.01) {
      return("**")
    } else if (p_value <= 0.05) {
      return("*")
    } else {
      return("")  # Not significant
    }
  }
  
  # Extract p-value from test results
  p_value = test_result$p.value
  
  # Convert p-value to stars
  stars = convert_p_to_stars(p_value)
  
  # Create a dataframe `annotations`
  annotations = data.frame(
    p_value = p_value,
    Stars = stars,
    stringsAsFactors = FALSE
  )
  
  # Filter out non-significant annotations
  annotations = annotations[annotations$Stars != "", ]
  
  return(annotations)
}

# Create annotations
annotation = prepare_annotations(test_result)

# Find the maximum y-values -> this is for dynamic annotation bars in the plots
annotation = annotation %>%
  mutate(
    max_y = max(group_means$Global_Mean + group_means$SD)
  )


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
plot_data = function(group_data, col_name, annotation_data, x = "Genotype", sd = "SD") {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(dapibtub_df[[col_name]], na.rm = TRUE)
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
    
    # Error bars
    geom_errorbar(aes(ymin = group_data[["Global_Mean"]] - group_data[[sd]],
                      ymax = group_data[["Global_Mean"]] + group_data[[sd]]),
                  width = 0.2,
                  position = position_dodge(0.9)) +
    
    # Graph titles
    labs(title = "",
         x = "",
         y = "β-III tubulin+ cells (%)",
         fill = x) +
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Define shapes for each DIFF value
  diff_shapes <- c("6" = 21, "7" = 22, "8" = 24)
  
  # Overlay individual data points with different shapes for DIFF
  p = p + geom_quasirandom(data = dapibtub_df, aes_string(x = x,
                                                              y = col_name,
                                                              shape = "DIFF"),  # Added shape aesthetic
                           width = 0.2, size = 1.5, fill = "black", alpha = 0.5) +  # Adjust width and alpha transparency here
    scale_shape_manual(values = diff_shapes)  # Adjust shape values if needed
  
  # Add conditional annotations for significant p-values
  if (nrow(annotation_data) > 0) {
    # Add significance stars
    # `x = 1.5` is used to position the text centrally between the two bars -> assumes that genotype has 2 ordered levels which correspond to position 1 and 2
    # `y = max_y * 1.08` places the text just above the estimated maximum y value
    # `position_dodge(width = 0.9` function is used to align the text with the corresponding bars
    p = p + geom_text(data = annotation_data, aes(label = Stars, x = 1.5, y = max_y_value * 1),
                      position = position_dodge(width = 0.9), inherit.aes = FALSE, vjust = -0.5,
                      size = 7)  # Adjust size here)
    
    # Add significance line
    # `x` and `xend` set the x-axis positions of the line
    # `y` and `yend` set the y-axis positions of the line
    # `position_dodge(width = 0.9` aligns the line with the bar positions
    p = p + geom_segment(data = annotation_data, aes(x = 1, xend = 2, y = max_y_value * 1.05, yend = max_y_value * 1.05),
                         linetype = "solid", color = "black", position = position_dodge(width = 0.9), inherit.aes = FALSE)
  }
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(group_means, "percentage_positive", annotation)

plot

# Export plot
# Open a PNG file to save the plot
#
# For plot title with 1 line:
# width=825, height=1335
#
# For plot title with 2 lines:
# width=825, height=1390
png(paste0(parent_filepath, "btub_percentage.png"), width=825, height=1335, res=300)

# Create a plot
plot

# Close the device
dev.off()

# Histograms
# plot2 = ggplot(dapibtub_df, aes(x = PRE_MeanIntensity, fill = DIFF, color = DIFF)) +
#   geom_density(alpha = 0.4, size = 0.5, adjust = 1.5) +  # Adjust the smoothness
#   labs(title = paste0("Density Plot of ", measurement, " (", marker, ")"),
#        x = measurement,
#        y = "Density") +
#   theme_minimal() +
#   theme(legend.position = "right") +
#   xlim(min(dapibtub_df$PRE_MeanIntensity) - 500, max(dapibtub_df$PRE_MeanIntensity) + 750) +
#   facet_wrap(~ Genotype, scales = "free_y")  # Separate plots by Genotype
# 
# plot2

# Export test results
# Function to extract p-value, method, alternative hypothesis, and sample sizes per group from test results
extract_test_results = function(test_results, data) {
  results_df = data.frame(
    p_value = numeric(),
    method = character(),
    alternative = character(),
    WT_sample_size = integer(),  # Sample size for WT genotype
    Q331K_sample_size = integer(),  # Sample size for Q331K genotype
    stringsAsFactors = FALSE
  )
  
  # Count the number of samples for each genotype within the current DIV
  wt_sample_count = nrow(subset(data, Genotype == "WT"))
  q331k_sample_count = nrow(subset(data, Genotype == "Q331K"))
  
  results_df = data.frame(
    p_value = round(test_result$p.value, 4),  # Round p-value to 4 decimal places
    method = test_result$method,
    alternative = test_result$alternative,
    WT_sample_size = wt_sample_count,  # Add WT sample size to the results
    Q331K_sample_size = q331k_sample_count  # Add Q331K sample size to the results
  )
  
  return(results_df)
}

# # Extract results
# results_df = extract_test_results(test_result, dapibtub_df)
# 
# # Define file path for saving CSVs
# csv_path = paste0(parent_filepath, relative_filepath, marker, "_puncta_", measurement, "intensity_raw.csv")
# 
# # Export the test results to CSV files
# write.csv(results_df, csv_path, row.names=FALSE)

# Print message
print(message)

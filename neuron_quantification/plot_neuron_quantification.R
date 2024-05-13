
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/echo/imaging_data_y1/"
relative_filepath = "neuron_quantification/"
filename = "neuron_quantification.xlsx"
marker = "MAP2"


# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
quantification_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Define cell_line variable as a factor with levels
quantification_data$Genotype = factor(quantification_data$Genotype, levels = c("WT", "Q331K"))


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

# Find mean for each sample
sample_means = mean_by_sample(quantification_data, "positive_cells/total_cells")


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


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_test = function(data) {
  
  # Extract values grouped by Genotype using N_Mean
  # Test for normality using Shapiro-Wilk test
  normality_result = by(data$N_Mean, data$Genotype, shapiro.test)
  
  wt_normal = normality_result$WT$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                 # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
  q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                       # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
  
  if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
    # If data is not normally distributed, perform Mann-Whitney test
    message = "Data not normally distributed; perform Mann-Whitney test"
    test_result = wilcox.test(N_Mean ~ Genotype, data = data) # `N_Mean ~ Genotype` specifies that you want to compare values in the `N_Mean` column grouped by `Genotype`
  } else {
    # Test for equal variances using F-test
    variance_test = var.test(N_Mean ~ Genotype, data = data)
    var_equal = variance_test$p.value > 0.05 # If p-value is greater than 0.05, it assumes both genotypes have equal variances
    
    if (!var_equal) { # If `var_equal` variable is set to `FALSE`
      # If variances are not equal, perform t-test with var.equal = FALSE
      message = "Data does not have equal variance; perform Welch's t-test"
      test_result = t.test(N_Mean ~ Genotype, data = data, var.equal = FALSE)
    } else {
      # If variances are equal, perform t-test with var.equal = TRUE
      message = "Data has equal variance; perform Student's t-test"
      test_result = t.test(N_Mean ~ Genotype, data = data, var.equal = TRUE)
    }
  }
  
  # Return both test results and the message
  return(list(test_result = test_result, message = message))
}

# Perform t-test
result = perform_test(sample_means)

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
          axis.title.y = element_text(margin = margin(r = 13), # Adjust y-axis title position
                                      size = 13), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, sample_data, annotation_data, x = "Genotype", sd = "SD") {
  # Check if the necessary columns exist in the group_data
  if (!all(c(x, y = "Global_Mean", sd) %in% names(group_data))) {
    stop("group_data does not contain the necessary columns.")
  }
  
  # Check if the necessary columns exist in the sample_data
  if (!all(c(x, y = "N_Mean") %in% names(sample_data))) {
    stop("sample_data does not contain the necessary columns.")
  }
  
  # Calculate the maximum y value to set upper axis limit
  # max_y_value = max(group_data["Global_Mean"] + group_data[sd], na.rm = TRUE)
  # upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
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
         y = paste0(marker, "+ cells (%)"),
         fill = x) +
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Overlay individual data points
  p = p + geom_point(data = sample_data, aes_string(x = x,
                                                    y = "N_Mean"),
                     position = position_dodge(0.9), size = 1.5)
  
  # Add conditional annotations for significant p-values
  if (nrow(annotation_data) > 0) {
    # Add significance stars
    # `x = 1.5` is used to position the text centrally between the two bars -> assumes that genotype has 2 ordered levels which correspond to position 1 and 2
    # `y = max_y * 1.08` places the text just above the estimated maximum y value
    # `position_dodge(width = 0.9` function is used to align the text with the corresponding bars
    p = p + geom_text(data = annotation_data, aes(label = Stars, x = 1.5, y = max_y * 1.08),
                      position = position_dodge(width = 0.9), inherit.aes = FALSE, vjust = -0.5,
                      size = 7)  # Adjust size here)
    
    # Add significance line
    # `x` and `xend` set the x-axis positions of the line
    # `y` and `yend` set the y-axis positions of the line
    # `position_dodge(width = 0.9` aligns the line with the bar positions
    p = p + geom_segment(data = annotation_data, aes(x = 1, xend = 2, y = max_y * 1.12, yend = max_y * 1.12),
                         linetype = "solid", color = "black", position = position_dodge(width = 0.9), inherit.aes = FALSE)
  }
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(group_means, sample_means, annotation)

# Export plot
# Open a PNG file to save the plot
#
# For plot title with 1 line:
# width=825, height=1335
#
# For plot title with 2 lines:
# width=825, height=1390
png(paste0(parent_filepath, relative_filepath, marker, "_quantification.png"), width=825, height=1335, res=300)

# Create a plot
plot

# Close the device
dev.off()

# Print message
print(message)

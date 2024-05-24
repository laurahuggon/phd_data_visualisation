
# For normalising to BTUB intensity: uncomment the calculation code, change the y-axis title on the plot, and change the file name.

# Load libraries ----------------------------------------------------------

library(tidyverse)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/cellbody_intensity/"

marker = "synaptophysin"
entity = "PRE"
measurement = "MedianIntensity"


# Load data ---------------------------------------------------------------

# Create relative_filepath and filename
if (marker == "synaptophysin") {
  relative_filepath = "syp/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 3), ".csv")
} else if (marker == "Homer") {
  relative_filepath = "hmr/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 3), ".csv")
} else if (marker == "syntaxin 1A") {
  relative_filepath = "stx/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 3), ".csv")
} else if (marker == "PSD-95") {
  relative_filepath = "psd/"
  filename = paste0(entity, "_cellbody_intensity_", substr(relative_filepath, 1, 3), ".csv")
}

full_filename = paste0(parent_filepath, relative_filepath, filename)
nis_elements_df = read_csv(full_filename)


# Prepare data ------------------------------------------------------------

# Unblind samples by adding Genotype, DIFF and DIV columns
# Define mappings - map each suffix to its corresponding value
if (marker == "synaptophysin" | marker == "PSD-95") {
  genotype_map = c("A" = "WT", "B" = "Q331K", "C" = "WT", "D" = "WT", "E" = "Q331K",
                   "F" = "Q331K", "G" = "Q331K", "H" = "WT", "I" = "Q331K", "J" = "WT",
                   "K" = "Q331K", "L" = "WT", "M" = "WT", "N" = "Q331K", "O" = "WT",
                   "P" = "Q331K", "Q" = "WT", "R" = "Q331K", "S" = "WT", "T" = "Q331K",
                   "U" = "WT", "V" = "Q331K", "W" = "Q331K", "X" = "WT")
  
  diff_map = c("A" = "3", "B" = "3", "C" = "3", "D" = "4", "E" = "3", "F" = "4", "G" = "4",
               "H" = "5", "I" = "3", "J" = "3", "K" = "5", "L" = "4", "M" = "4", "N" = "3",
               "O" = "3", "P" = "4", "Q" = "5", "R" = "4", "S" = "5", "T" = "5", "U" = "4",
               "V" = "5", "W" = "5", "X" = "5")
  
  div_map = c("A" = "14", "B" = "7", "C" = "7", "D" = "7", "E" = "14", "F" = "7",
              "G" = "14", "H" = "7", "I" = "21", "J" = "21", "K" = "7", "L" = "14",
              "M" = "21", "N" = "28", "O" = "28", "P" = "21", "Q" = "21", "R" = "28",
              "S" = "14", "T" = "21", "U" = "28", "V" = "14", "W" = "28", "X" = "28")
  
  suffix_regex = "(?<=_)[A-X](?=_)" # Uses lookbehind `(?<=_)` and lookahead `(?=_)` to capture character between two underscores
  
} else if (marker == "syntaxin 1A" | marker == "Homer") {
  genotype_map = c("1" = "Q331K", "2" = "Q331K", "3" = "WT", "4" = "WT",
                   "5" = "Q331K", "6" = "WT", "7" = "WT", "8" = "Q331K",
                   "9" = "Q331K", "10" = "WT", "11" = "Q331K", "12" = "Q331K",
                   "13" = "WT", "14" = "Q331K", "15" = "WT", "16" = "WT",
                   "17" = "WT", "18" = "Q331K", "19" = "Q331K", "20" = "WT",
                   "21" = "WT", "22" = "Q331K", "23" = "WT", "24" = "Q331K",
                   "25" = "Q331K", "26" = "WT", "27" = "Q331K", "28" = "WT",
                   "29" = "WT", "30" = "Q331K", "31" = "WT", "32" = "Q331K")
  
  diff_map = c("1" = "4", "2" = "3", "3" = "5", "4" = "3", "5" = "4",
               "6" = "4", "7" = "3", "8" = "5", "9" = "3", "10" = "4",
               "11" = "5", "12" = "3", "13" = "5", "14" = "3", "15" = "3",
               "16" = "3", "17" = "4", "18" = "4", "19" = "5", "20" = "4",
               "21" = "5", "22" = "4", "23" = "5", "24" = "5", "25" = "14",
               "26" = "14", "27" = "14", "28" = "14", "29" = "14", "30" = "14",
               "31" = "14", "32" = "14")
  
  div_map = c("1" = "7", "2" = "14", "3" = "7", "4" = "7", "5" = "14",
              "6" = "7", "7" = "14", "8" = "7", "9" = "7", "10" = "14",
              "11" = "14", "12" = "21", "13" = "14", "14" = "28", "15" = "21",
              "16" = "28", "17" = "21", "18" = "21", "19" = "21", "20" = "28",
              "21" = "21", "22" = "28", "23" = "28", "24" = "28", "25" = "7",
              "26" = "7", "27" = "14", "28" = "14", "29" = "21", "30" = "21",
              "31" = "28", "32" = "28")
  
  suffix_regex = "(?<=_)[0-9]+(?=_[0-9]{3})" # Matches a sequence of one or more digits `[0-9]+` that are preceded by an underscore `(?<=_)` and followed by an underscore and exactly three digits `(?=_[0-9]{3})`
  
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

# Define Genotype and DIV variable as a factor with levels
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
          axis.title.y = element_text(margin = margin(r = 7.5), # Adjust y-axis title position
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
  max_y_value = max(group_data["Global_Mean"] + group_data[sd], na.rm = TRUE)
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

# Print message
print(message)

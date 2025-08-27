
# For normalising intensities to BTUB intensity: uncomment the calculation code, add "Normalised" the y-axis title on the plot, and change the file name from "_raw" to "_normalised".

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggbeeswarm)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/global_intensity/"

marker = "Homer-1"
normalise = TRUE

# Load data ---------------------------------------------------------------

# Create relative_filepath and filename
if (marker == "synaptophysin") {
  filename = "syp/PRE_global_intensity_syp.csv"
} else if (marker == "Homer-1") {
  filename = "hmr/POST_global_intensity_hmr.csv"
} else if (marker == "syntaxin-1A") {
  filename = "stx/PRE_global_intensity_stx.csv"
} else if (marker == "PSD-95") {
  filename = "psd/POST_global_intensity_psd.csv"
}

full_filename = paste0(parent_filepath, filename)
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
  
} else if (marker == "syntaxin-1A" | marker == "Homer-1") {
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

# Remove DIFF 5
if(marker == "syntaxin-1A" | marker == "Homer-1") {
  nis_elements_df = nis_elements_df %>%
    filter(DIFF != "5")
}

# Define Genotype and DIFF variable as a factor with levels
nis_elements_df$Genotype = factor(nis_elements_df$Genotype, levels = c("WT", "Q331K"))
nis_elements_df$DIFF = factor(nis_elements_df$DIFF, levels = c(3, 4, 5, 14))

# Normalise MeanIntensity to BTUB_MeanIntensity
if(normalise == TRUE && (marker == "synaptophysin" | marker == "syntaxin-1A")) {
  nis_elements_df = nis_elements_df %>%
    mutate(
      mean_intensity = PRE_MeanIntensity / BTUB_MeanIntensity
    )
} else if (normalise == TRUE && (marker == "PSD-95" | marker == "Homer-1")) {
  nis_elements_df = nis_elements_df %>%
    mutate(
      mean_intensity = POST_MeanIntensity / BTUB_MeanIntensity
    )
} else if (normalise == FALSE && (marker == "synaptophysin" | marker == "syntaxin-1A")) {
  nis_elements_df = nis_elements_df %>%
    mutate(
      mean_intensity = PRE_MeanIntensity
    )
} else {
  nis_elements_df = nis_elements_df %>%
    mutate(
      mean_intensity = POST_MeanIntensity
    )
}

# Find sample means -------------------------------------------------------

# Create function that finds sample means for a given variable
mean_by_sample = function(data, column_name) {
  # Group by Genotype, DIFF
  result = data %>%
    group_by(Genotype, DIFF) %>%
    summarise(
      N = n(), # Count the number of images in each sample
      N_Mean = mean(.data[[column_name]], na.rm=TRUE), # Calculate mean for each sample
    )
  return(result)
}

# Find mean for each sample
sample_means = mean_by_sample(nis_elements_df, "mean_intensity")


# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, column_name) {
  # Group by Genotype
  result = data %>%
    group_by(Genotype) %>%
    summarise(
      N = n(), # Count the number of samples in each group
      Group_Mean = mean(.data[[column_name]]), # Calculate mean for each group
      SD = sd(.data[[column_name]]) # Calculate the SD for each group
    )
  return(result)
}

# Find mean for each group
group_means = mean_by_group(sample_means, "N_Mean")


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_test = function(data, value_col, grouping_col, ctrl_group, exp_group) {
  
  # Extract data for groups
  # [[1]] pulls the first (and only) column as a vector
  ctrl_data = data[data[[grouping_col]] == ctrl_group, value_col][[1]]
  exp_data = data[data[[grouping_col]] == exp_group, value_col][[1]]
  
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
test_result_list = perform_test(sample_means, "N_Mean", "Genotype", "WT", "Q331K")

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
  mutate(max_y = max(nis_elements_df$mean_intensity)) %>%
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
          strip.text = element_text(size = 12, face = "bold") # Facet title size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, individual_data, individual_col_name, x, test_results = test_result) {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(individual_data[[individual_col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Create y-axis title
  if(normalise == TRUE) {
    y_title = paste0("Normalised mean intensity\nof ", marker)
  } else {
    y_title = paste0("Mean intensity\nof ", marker, " (a.u.)")
  }
  
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
    labs(title = "Global",
         x = "",
         y = y_title,
         fill = x) + # Legend title
    
    # Plot appearance
    my_theme()
  
  if(normalise == TRUE) {
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0), labels = function(x) sprintf("%.2f", x))  # Setting both multiplier and add-on to 0
  } else {
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  }
  
  # Define shapes for each DIFF value
  diff_shapes = c("3" = 21, "4" = 22, "5" = 24, "14" = 25)
  
  # Overlay individual data points (optional with different shapes for DIFF)
  p = p + geom_point(data=individual_data, aes_string(x = x,
                                                      y = individual_col_name,
                                                      shape = "DIFF"),
                     size = 1.25,
                     fill = "black") +
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
plot = plot_data(group_means, "Group_Mean", sample_means, "N_Mean", "Genotype")

plot

if(normalise == TRUE) {
  # Define directory
  dir = sub("^([^/]+/).*", "\\1", filename)
  # Save plot
  ggsave(paste0(parent_filepath, dir, marker, "_global_mean-intensity_norm_pooled.png"), plot=plot, width=2.05, height=3.5, dpi=300, bg="white")
  # Export test results
  # Define file path for saving CSVs
  csv_path = paste0(parent_filepath, dir, "stats/", marker, "_global_mean-intensity_norm_pooled.csv")
  # Export the test results to CSV files
  write.csv(test_result, csv_path, row.names=FALSE)
} else {
  # Define directory
  dir = sub("^([^/]+/).*", "\\1", filename)
  # Save plot
  ggsave(paste0(parent_filepath, dir, marker, "_global_mean-intensity_raw_pooled.png"), plot=plot, width=2.1, height=3.5, dpi=300, bg="white")
  # Export test results
  # Define file path for saving CSVs
  csv_path = paste0(parent_filepath, dir, "stats/", marker, "_global_mean-intensity_raw_pooled.csv")
  # Export the test results to CSV files
  write.csv(test_result, csv_path, row.names=FALSE)
}


# Grouping of samples -----------------------------------------------------

# My data: each of my 24 samples (e.g. A, B, C, ...) has 10 images (e.g. 001, 002, 003, ...)
# NIS-Elements returns the mean values for each image (e.g. mean of A_001) = 240 values
# mean_by_sample returns the mean of the 10 images for each sample (e.g. mean of A) = 24 values
# mean_by_group returns the mean of the 3 samples in each group (defined by genotype and age) (e.g. mean of A, B and C) = 8 values


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(scales)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/synapse_morphology/"
relative_filepath = "n_1-3/syp_coloc/"
filename = "COLOCPRE_BTUB_A-X.csv"
entity = "PRE"
pre_marker = "synaptophysin"
post_marker = "PSD-95"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
nis_elements_df = read_csv(full_filename)


# Prepare data ------------------------------------------------------------

# Unblind samples by adding Genotype, DIFF and DIV columns
# Define mappings - map each suffix to its corresponding value
if (pre_marker == "synaptophysin" && post_marker == "PSD-95") {
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
  
} else if (pre_marker == "syntaxin-1A" && post_marker == "Homer-1") {
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
  mutate(Suffix = str_extract(Filename, suffix_regex)) %>% # Uses lookbehind `(?<=_)` and lookahead `(?=_)` to capture character between two underscores
  mutate(
    Genotype = genotype_map[Suffix], # Map suffixes to respective values using predefined mappings
    DIFF = diff_map[Suffix],
    DIV = div_map[Suffix]
  ) %>%
  select(-Suffix)  # Remove the Suffix column if it's not needed later

# Remove DIV 28
nis_elements_df = nis_elements_df %>%
  filter(DIV != "28")

# Remove DIFF 5
if(startsWith(relative_filepath, "n_1_2_4")) {
  nis_elements_df = nis_elements_df %>%
    filter(DIFF != "5")
}

# Define Genotype and DIV variable as a factor with levels
nis_elements_df$Genotype = factor(nis_elements_df$Genotype, levels = c("WT", "Q331K"))
nis_elements_df$DIV = factor(nis_elements_df$DIV, levels = c("7", "14", "21"))


# Find sample means -------------------------------------------------------

# Create function that finds sample means for a given variable
# mean_by_sample = function(data, column_name) {
#   # Group by Genotype, DIV, DIFF
#   result = data %>%
#     group_by(Genotype, DIV, DIFF) %>%
#     summarise(
#       N = n(), # Count the number of images in each sample
#       N_Mean = mean(.data[[column_name]]), # Calculate mean for each sample
#     )
#   return(result)
# }
# 
# # Find mean colocalisation, density and volume for each sample
# sample_mean_coloc = mean_by_sample(nis_elements_df, "Coloc")
# sample_mean_density = mean_by_sample(nis_elements_df, "DensityColoc")
# sample_mean_volume = mean_by_sample(nis_elements_df, "MeanVolumeColoc")

# Find group means --------------------------------------------------------

# Create function that finds group means for a given variable
mean_by_group = function(data, column_name) {
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
group_mean_coloc = mean_by_group(nis_elements_df, "Coloc")
group_mean_density = mean_by_group(nis_elements_df, "DensityColoc")
group_mean_volume = mean_by_group(nis_elements_df, "MeanVolumeColoc")


# Create a function that calculates the maximum y-values per facet -> this is for dynamic annotation bars in the plots
calculate_max_y_per_facet = function(data=nis_elements_df, col_name) {
  # Add a new column to the group_data that calculates the potential max height for the annotation bars
  data$max_y = data[[col_name]]
  
  # Aggregate these max heights by DIV to get the maximum for each facet
  max_y_per_div = aggregate(max_y ~ DIV, data = data, max)
  
  return(max_y_per_div)
}

# Find maximum y-values per facet for each data type
max_y_per_div_coloc = calculate_max_y_per_facet(col_name="Coloc")
max_y_per_div_density = calculate_max_y_per_facet(col_name="DensityColoc")
max_y_per_div_volume = calculate_max_y_per_facet(col_name="MeanVolumeColoc")


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_tests_with_report = function(data, col_name) {
  # Extract unique DIV levels to perform the t-test for each level separately
  div_levels = unique(data$DIV) # Extract unique values from the DIV column
  
  # Create an empty list to store test results
  test_results = list()
  
  # Initialize the report with the 'Message' column
  report = data.frame(Message = character(), stringsAsFactors = FALSE)
  
  # Loop through each DIV level - allows performance of separate analyses at each time point
  for (div in div_levels) {
    # Subset data by DIV
    data_div = subset(data, DIV == div) # Extract rows from dataframe where DIV matches the current DIV level in the loop
    
    # Extract Mean values for each genotype
    wt_data = subset(data_div, Genotype == "WT")[[col_name]]
    q331k_data = subset(data_div, Genotype == "Q331K")[[col_name]]
    
    # Test for normality using Shapiro-Wilk test
    wt_normal = shapiro.test(wt_data)$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
    # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
    q331k_normal = shapiro.test(q331k_data)$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
    # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
    
    if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
      # If either group is not normally distributed, perform Mann-Whitney test
      message = paste("DIV", div, ": data not normally distributed; perform Mann-Whitney test")
      report = rbind(report, data.frame(Message = message, stringsAsFactors = FALSE))
      test_results[[div]] = wilcox.test(wt_data, q331k_data) # `wt_data` and `q331k_data` are vectors containing the `Mean` values for each genotype at the current DIV level
      # Store the results in a list, `test_results`, under a key created using the current DIV value (`div`)
    } else {
      # If both groups are normally distributed, test for equal variance using the F-test
      var_equal = var.test(wt_data, q331k_data)$p.value > 0.05 # If p-value is greater than 0.05, it assumes both genotypes have equal variances
      
      if (!var_equal) { # If `var_equal` variable is set to `FALSE`
        # If variances are not equal, perform t-test with var.equal = FALSE
        message = paste("DIV", div, ": data does not have equal variance; perform Welch's t-test")
        report = rbind(report, data.frame(Message = message, stringsAsFactors = FALSE))
        test_results[[div]] = t.test(wt_data, q331k_data, var.equal = FALSE) 
      } else {
        # If variances are equal, perform t-test with var.equal = TRUE
        test_results[[div]] = t.test(wt_data, q331k_data, var.equal = TRUE)
      }
    }
  }
  
  # Return both test results and the report
  return(list(test_results = test_results, report = report))
}

# Perform the tests with reporting for each data type
coloc_results = perform_tests_with_report(nis_elements_df, "Coloc")
density_results = perform_tests_with_report(nis_elements_df, "DensityColoc")
volume_results = perform_tests_with_report(nis_elements_df, "MeanVolumeColoc")

# Accessing the test results
coloc_test_results = coloc_results$test_results
density_test_results = density_results$test_results
volume_test_results = volume_results$test_results

# Collect all the reports into a named list
reports_list = list(
  Coloc = coloc_results$report,
  Density = density_results$report,
  Volume = volume_results$report
)

# Create function that converts p-values to signficance stars and generates a dataframe of annotations
prepare_annotations = function(test_results) {
  # Convert p-values to stars based on traditional significance levels
  convert_p_to_stars = function(p_value) {
    if (p_value < 0.0001) {
      return("****")
    } else if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else {
      return("")  # Not significant
    }
  }
  
  # Create a dataframe `annotations`
  annotations = data.frame(
    DIV = names(test_results), # One entry per DIV level
    p_value = sapply(test_results, function(x) x$p.value), # Extract p-values
    stringsAsFactors = FALSE
  )
  
  # Convert p-values to stars
  annotations$Stars = sapply(annotations$p_value, convert_p_to_stars)
  
  # Set factor levels for `DIV`
  annotations$DIV = factor(annotations$DIV, levels = c("7", "14", "21"))
  
  # Filter out non-significant annotations
  annotations = annotations[annotations$Stars != "", ]
  
  return(annotations)
}

# Create annotations for each data type
coloc_annotations = prepare_annotations(coloc_test_results)
density_annotations = prepare_annotations(density_test_results)
volume_annotations = prepare_annotations(volume_test_results)

# Create a function that merges max y-values with annotations -> this is for dynamic annotation bars in the plots
merge_annotations_with_max_y = function(annotations, max_y_per_div) {
  # Merge the max y-values per DIV with the annotations data frame
  annotations = merge(annotations, max_y_per_div, by = "DIV")
  
  return(annotations)
}

# Merge max y-values with annotations for each data type
coloc_annotations = merge_annotations_with_max_y(coloc_annotations, max_y_per_div_coloc)
density_annotations = merge_annotations_with_max_y(density_annotations, max_y_per_div_density)
volume_annotations = merge_annotations_with_max_y(volume_annotations, max_y_per_div_volume)


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
          axis.title.y = element_text(margin = margin(r = 15), # Adjust y-axis title position
                                      size = 12), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10), # Increase y-axis text size
          # Facet-specific
          panel.spacing = unit(1, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 11, face = "bold")  # Increase facet title size and make it bold
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_by_genotype_div = function(group_data, col_name, annotation_data, x = "Genotype", sd = "SD", facet = "DIV") {
  # Extracting the name of the dataframe
  data_name = deparse(substitute(group_data)) # Get the name of the `group_data` dataframe as a string
  
  # Determine y-axis label based on the inferred dataframe name
  if (grepl("coloc", data_name)) {
    y_label = paste0("Colocalised ", protein_name, "\npuncta (%)")
  } else if (grepl("density", data_name)) {
    y_label = paste0("Density of ", protein_name, "\npuncta (puncta/μm³)")
  } else if (grepl("volume", data_name)) {
    y_label = paste0("Volume of ", protein_name, "\npuncta (μm³)")
  } else {
    y_label = "Mean"
  }
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(nis_elements_df[[col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Define custom labeller function to add "DIV " to facet titles
  my_labeller = as_labeller(function(x) paste("DIV", x))
  
  # Create the bar plot
  p = ggplot(group_data, aes_string(x = x,
                                    y = "Global_Mean",
                                    fill = x)) +
    # Bar plot
    geom_col(position = position_dodge(0.9),
             width = 0.8,
             color = "black") +
    scale_fill_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF")) +
    
    # Error bars
    geom_errorbar(aes(ymin = group_data[["Global_Mean"]] - group_data[[sd]],
                      ymax = group_data[["Global_Mean"]] + group_data[[sd]]),
                  width = 0.2,
                  position = position_dodge(0.9)) +
    
    # Facet
    facet_wrap(as.formula(paste0("~ ", facet)), ncol = 3, labeller = my_labeller, axes="all") +
    
    # Graph titles
    labs(x = "",
         y = y_label,
         fill = x) +
    
    # Plot appearance
    my_theme_facet()
  
  # Set y-axis ticks
  if (grepl("coloc", data_name)) {
    upper_limit = 110
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))
  } else {
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0), labels = label_number(accuracy=0.01))
  }
  # Define shapes for each DIFF value
  diff_shapes <- c("3" = 21, "4" = 22, "5" = 24, "14" = 25)
  
  # Overlay individual data points with different shapes for DIFF
  p = p + geom_quasirandom(data = nis_elements_df, aes_string(x = x,
                                                              y = col_name,
                                                              shape = "DIFF"),  # Added shape aesthetic
                           width = 0.2, size = 1.25, fill = "black", alpha = 0.5) +  # Adjust width and alpha transparency here
    scale_shape_manual(values = diff_shapes)  # Adjust shape values if needed
  
  # Add conditional annotations for significant p-values
  # The placement of annotations specific to facets relies on the use of `annotation_data` that contains significance annotations and max y-values that are merged with specific facet information (`DIV`)
  # Uses the `DIV` column in `annotation_data` to apply significant annnotations to the corresponding facets
  if (nrow(annotation_data) > 0) {
    # Add significance stars
    # `x = 1.5` is used to position the text centrally between the two bars -> assumes that genotype has 2 ordered levels which correspond to position 1 and 2
    # `y = max_y * 1.08` places the text just above the estimated maximum y value
    # `position_dodge(width = 0.9` function is used to align the text with the corresponding bars
    p = p + geom_text(data = annotation_data, aes(label = Stars, x = 1.5, y = max_y + 0.05*upper_limit),
                      position = position_dodge(width = 0.9), inherit.aes = FALSE, vjust = -0.5,
                      size = 6)  # Adjust size here
    
    # Add significance line
    # `x` and `xend` set the x-axis positions of the line
    # `y` and `yend` set the y-axis positions of the line
    # `position_dodge(width = 0.9` aligns the line with the bar positions
    p = p + geom_segment(data = annotation_data, aes(x = 1, xend = 2, y = (max_y + 0.09*upper_limit), yend = (max_y + 0.09*upper_limit)),
                         linetype = "solid", color = "black", position = position_dodge(width = 0.9), inherit.aes = FALSE)
  }
  
  # Print the plot
  return(p)
}

# Make plots for each data type
coloc_plot = plot_by_genotype_div(group_mean_coloc, "Coloc", coloc_annotations)
density_plot = plot_by_genotype_div(group_mean_density, "DensityColoc", density_annotations)
volume_plot = plot_by_genotype_div(group_mean_volume, "MeanVolumeColoc", volume_annotations)

# Arrange plots
all_plots = coloc_plot /
  density_plot /
  volume_plot

all_plots

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_coloc_all_plots_allimages.png"), plot=all_plots, width=5.25, height=10.5, dpi=300, bg="white")

# Export test results
# Function to extract p-value, method, alternative hypothesis, and sample sizes per group from test results
extract_test_results = function(test_results, data) {
  results_df = data.frame(
    DIV = character(),
    p_value = numeric(),
    method = character(),
    alternative = character(),
    WT_sample_size = integer(),  # Sample size for WT genotype
    Q331K_sample_size = integer(),  # Sample size for Q331K genotype
    stringsAsFactors = FALSE
  )
  
  for (div in names(test_results)) {
    test = test_results[[div]]
    
    # Count the number of samples for each genotype within the current DIV
    wt_sample_count = nrow(subset(data, DIV == div & Genotype == "WT"))
    q331k_sample_count = nrow(subset(data, DIV == div & Genotype == "Q331K"))
    
    new_row = data.frame(
      DIV = div,
      p_value = round(test$p.value, 4),  # Round p-value to 4 decimal places
      method = test$method,
      alternative = test$alternative,
      WT_sample_size = wt_sample_count,  # Add WT sample size to the results
      Q331K_sample_size = q331k_sample_count  # Add Q331K sample size to the results
    )
    
    results_df = rbind(results_df, new_row)
  }
  
  return(results_df)
}

# Extract results for colocalization, density, and volume
coloc_test_results_df = extract_test_results(coloc_test_results, nis_elements_df)
density_test_results_df = extract_test_results(density_test_results, nis_elements_df)
volume_test_results_df = extract_test_results(volume_test_results, nis_elements_df)

# Convert DIV to a factor with specified levels
coloc_test_results_df$DIV <- factor(coloc_test_results_df$DIV, levels = c("7", "14", "21"), ordered = TRUE)
density_test_results_df$DIV <- factor(density_test_results_df$DIV, levels = c("7", "14", "21"), ordered = TRUE)
volume_test_results_df$DIV <- factor(volume_test_results_df$DIV, levels = c("7", "14", "21"), ordered = TRUE)

# Order by DIV
coloc_test_results_df <- coloc_test_results_df[order(coloc_test_results_df$DIV), ]
density_test_results_df <- density_test_results_df[order(density_test_results_df$DIV), ]
volume_test_results_df <- volume_test_results_df[order(volume_test_results_df$DIV), ]

# Define file paths for saving CSVs
coloc_csv_path = paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_coloc_test_results_allimages.csv")
density_csv_path = paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_density_test_results_allimages.csv")
volume_csv_path = paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_volume_test_results_allimages.csv")

# Export the test results to CSV files
write.csv(coloc_test_results_df, coloc_csv_path, row.names = FALSE)
write.csv(density_test_results_df, density_csv_path, row.names = FALSE)
write.csv(volume_test_results_df, volume_csv_path, row.names = FALSE)

# Loop through the list and print only non-empty reports
for (name in names(reports_list)) {
  report = reports_list[[name]]
  if (nrow(report) > 0) {
    # Add a title to distinguish reports
    cat(paste("\n", name, "Report:\n"))
    # Print only the message content without headers
    cat(report$Message, sep = "\n")
  }
}


# Density plot ------------------------------------------------------------

my_labeller = as_labeller(function(x) paste("DIV", x))

kde = ggplot(data=nis_elements_df, aes(x=MeanVolumeColoc, color=Genotype, fill=Genotype)) +
  geom_density(alpha=0.2) +
  facet_wrap(~ DIV, labeller = my_labeller, ncol=3, axes="all") +
  scale_color_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF")) +
  scale_fill_manual(values = c("WT" = "#F3D99E", "Q331K" = "#DBAEAF")) +
  geom_vline(data=group_mean_volume, aes(xintercept=Global_Mean, colour=Genotype), linetype="dashed") +
  my_theme_facet() +
  theme(legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 10), size = 12)) + # Adjust x-axis title position
  labs(x = paste0("Mean volume of ", protein_name, " puncta (μm³)"))

kde

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_coloc_density.png"), plot=kde, width=14, height=3, dpi=300, bg="white")
                   
                  
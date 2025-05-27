
# Grouping of samples -----------------------------------------------------

# My data: each of my 24 samples (e.g. A, B, C, ...) has 10 images (e.g. 001, 002, 003, ...)
# NIS-Elements returns the mean values for each image (e.g. mean of A_001) = 240 values
# mean_by_sample returns the mean of the 10 images for each sample (e.g. mean of A) = 24 values
# mean_by_group returns the mean of the 3 samples in each group (defined by genotype and age) (e.g. mean of A, B and C) = 8 values


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/synapse_morphology/"
relative_filepath = "n_1_2_4/stx_coloc/"
filename = "COLOCPRE_BTUB_1-32.csv"
entity = "PRE"
pre_marker = "syntaxin-1A"
post_marker = "Homer-1"


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

# Split dataframe per metric
coloc = nis_elements_df %>%
  select(Filename:Entity, ColocCount:Coloc, Genotype:DIV)

density = nis_elements_df %>%
  select(Filename:Entity, ColocCount, TotalNeuriteVolume, DensityColoc, Genotype:DIV)

volume = nis_elements_df %>%
  select(Filename:Entity, MeanVolumeColoc:ColocCount, Genotype:DIV)


# Find sample means -------------------------------------------------------

mean_by_sample = function(data, value_col, grouping_col) {
  # Group and summarise
  result = data %>%
    group_by(!!!syms(grouping_col)) %>%
    summarise(
      N = n(), # Count the number of images in each sample
      Sample_Mean = mean(.data[[value_col]]), # Calculate mean for each sample
    )
  return(result)
}

# Find mean for each sample
sample_mean_coloc = mean_by_sample(coloc, "Coloc", c("DIV", "Genotype", "DIFF"))
sample_mean_density = mean_by_sample(density, "DensityColoc", c("DIV", "Genotype", "DIFF"))
sample_mean_volume= mean_by_sample(volume, "MeanVolumeColoc", c("DIV", "Genotype", "DIFF"))

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
group_mean_coloc = mean_by_group(sample_mean_coloc, "Sample_Mean", c("DIV", "Genotype"))
group_mean_density = mean_by_group(sample_mean_density, "Sample_Mean", c("DIV", "Genotype"))
group_mean_volume = mean_by_group(sample_mean_volume, "Sample_Mean", c("DIV", "Genotype"))


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

create_test_result_df = function(data, facet_grouping) {
  # group_by and summarise can be used to iterate over groups of data
  test_result = data %>%
    # Group by facet_grouping column, which creates a subset of long_data, where each subset corresponds to a facet group
    group_by(!!!syms(facet_grouping)) %>%
    # Apply the perform_test function to each subset of data
    summarise(
      n_ctrl = do.call(perform_test, c(list(data=cur_data()), args))$n_ctrl,
      n_exp = do.call(perform_test, c(list(data=cur_data()), args))$n_exp,
      p.value = do.call(perform_test, c(list(data=cur_data()), args))$test_result$p.value,
      method = do.call(perform_test, c(list(data=cur_data()), args))$test_result$method,
      alternative = do.call(perform_test, c(list(data=cur_data()), args))$test_result$alternative
    )
  
  # Apply multiple comparisons adjustment
  # test_result$p.adjust = p.adjust(test_result$p.value, method = "bonferroni")
  
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
  
  # Find maximum y-value per facet
  max_y_per_facet = data %>%
    group_by(!!!syms(facet_grouping)) %>%
    summarise(
      max_y = max(.data[[args$value_col]])
    )
  
  # Merge maximum y-values with test_result
  test_result = merge(test_result, max_y_per_facet, by=facet_grouping)
  
  # Replace max_y with NA if not significant
  test_result = test_result %>%
    mutate(max_y = ifelse(stars == "", NA, max_y)) %>%
    arrange(!!!syms(facet_grouping))
  
  return(test_result)
  
}

# Perform t-test for each facet
args = list(
  value_col = "Sample_Mean",
  grouping_col = "Genotype", 
  ctrl_group = "WT",
  exp_group = "Q331K"
)

test_result_coloc = create_test_result_df(sample_mean_coloc, "DIV")
test_result_density = create_test_result_df(sample_mean_density, "DIV")
test_result_volume = create_test_result_df(sample_mean_volume, "DIV")

# Data visualisation ------------------------------------------------------

# Define the protein name for y-axis label
protein_name = if (entity == "PRE") {
  pre_marker
} else {
  post_marker
}

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold", hjust = 0.5), # Adjust plot title
          axis.title.y = element_text(margin = margin(r = 15), size = 12), # Adjust y-axis title
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10), # Increase y-axis text size
          # Facet-specific
          panel.spacing = unit(0.5, "lines"), # Adjust spacing between facet panels
          strip.text = element_text(size = 10, face="bold") # Facet title size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, individual_data, individual_col_name, x, facet_grouping, test_results = test_result) {
  
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
  max_y_value = max(individual_data[[individual_col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Define custom labeller function to add "DIV " to facet titles
  my_labeller = as_labeller(function(x) paste("DIV", x))
  
  # Reformulate
  facet_formula = reformulate(facet_grouping)
  
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
    
    # Facet
    facet_wrap(facet_formula, nrow=1, labeller = my_labeller, axes= "all") +
    
    # Graph titles
    labs(title = "",
         x = "",
         y = y_label,
         fill = x) + # Legend title
    
    # Plot appearance
    my_theme()
  
  if (grepl("coloc", data_name)) {
    p = p + scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + # Setting both multiplier and add-on to 0
      theme(axis.title.y = element_text(margin = margin(r = 17.75)),
        panel.spacing = unit(0.7, "lines"))
  } else {
    p = p + scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0), labels = label_number(accuracy=0.01))  # Setting both multiplier and add-on to 0
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
                    position = position_dodge(width = 0.5),
                    inherit.aes = FALSE,
                    size = 6)  # Adjust size here
  # Significance lines
  p = p + geom_segment(data = test_results, aes(x = 1,
                                                xend = 2,
                                                y = (max_y + 0.075*upper_limit),
                                                yend = (max_y + 0.075*upper_limit)),
                       linetype = "solid",
                       color = "black",
                       position = position_dodge(width = 0.5),
                       inherit.aes = FALSE)
  
  # Print the plot
  return(p)
}

# Make plot
coloc_plot = plot_data(group_mean_coloc, "Group_Mean", sample_mean_coloc, "Sample_Mean", "Genotype", "DIV", test_results = test_result_coloc)
density_plot = plot_data(group_mean_density, "Group_Mean", sample_mean_density, "Sample_Mean", "Genotype", "DIV", test_results = test_result_density)
volume_plot = plot_data(group_mean_volume, "Group_Mean", sample_mean_volume, "Sample_Mean", "Genotype", "DIV", test_results = test_result_volume)
volume_plot

# Save plot
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_coloccoloc.png"), plot=coloc_plot, width=4.65, height=3.5, dpi=300, bg="white")
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_colocdensity.png"), plot=density_plot, width=4.65, height=3.5, dpi=300, bg="white")
ggsave(paste0(parent_filepath, relative_filepath, protein_name, "_colocvolume.png"), plot=volume_plot, width=4.65, height=3.5, dpi=300, bg="white")

# Export the test results to CSV files
write.csv(test_result_coloc, paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_test_result_coloc.csv"), row.names=FALSE)
write.csv(test_result_density, paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_test_result_density.csv"), row.names=FALSE)
write.csv(test_result_volume, paste0(parent_filepath, relative_filepath, "stats/", protein_name, "_test_result_volume.csv"), row.names=FALSE)


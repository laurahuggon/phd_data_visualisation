
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/mouse_cytosol/"
filename = "cytosol_allproteins.xlsx"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Define Genotype as a factor with levels
empiria_data$Genotype = factor(empiria_data$Genotype, levels = c("NTg", "Q331K"))

# Reshape data into long format
long_data = empiria_data %>%
  select(-Name) %>%
  pivot_longer(-c(Genotype, Sex), names_to = "Protein", values_to = "Signal") %>%
  arrange(Genotype) %>%
  arrange(Protein)

# Define Protein as a factor with levels
long_data$Protein = factor(long_data$Protein, levels = c("Synaptophysin", "Syntaxin-1A", "VAMP2", "SNAP-25", "Munc13-1", "Munc18-1", "PSD-95", "Homer-1"))


# Normalise to control mean -----------------------------------------------

# Use long_data to calculate to calculate mean for each genotype for each protein
group_means = long_data %>%
  group_by(Genotype, Protein) %>%
  summarise(
    mean = mean(Signal),
    .groups = "drop"
  )

# Extract the mean of the control group for each protein
control_means = group_means %>%
  filter(Genotype == "NTg") %>%
  select(Protein, control_mean = mean)

# Divide each sample by the control mean (per protein)
long_data = long_data %>%
  left_join(control_means, by = "Protein") %>%
  mutate(norm_signal = Signal / control_mean) %>%
  select(Genotype:Signal, norm_signal) %>%
  mutate(norm_signal = ifelse(is.na(norm_signal), 0, norm_signal))

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
group_means = mean_by_group(long_data, "norm_signal", c("Protein", "Genotype"))


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
  value_col = "norm_signal",
  grouping_col = "Genotype", 
  ctrl_group = "NTg",
  exp_group = "Q331K"
)

# test_result = create_test_result_df(long_data, "Protein")
test_result = create_test_result_df(long_data %>% filter(!Protein %in% c("Munc13-1", "Munc18-1")), "Protein")


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(
      axis.line = element_line(colour = "black"),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12), # Adjust plot title
      axis.title.x = element_text(margin = margin(t = 15), size = 12), # Adjust x-axis title
      axis.title.y = element_text(margin = margin(r = 15), size = 12), # Adjust y-axis title
      axis.text.x = element_text(size = 10), # Increase x-axis text size
      axis.text.y = element_text(size = 10), # Increase y-axis text size
      # Facet-specific
      panel.spacing = unit(0.5, "lines"), # Adjust spacing between facet panels
      strip.text = element_text(size = 10) # Facet title size
    ) 
}


# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_data = function(group_data, group_col_name, individual_data, individual_col_name, x, facet_grouping, test_results = test_result) {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(individual_data[[individual_col_name]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 25% buffer above the max value
  
  # Reformulate
  facet_formula = reformulate(facet_grouping)
  
  # Create the bar plot
  p = ggplot(group_data, aes(x = .data[[x]], 
                             y = .data[[group_col_name]],
                             fill = .data[[x]])) +
    
    # Bar plot
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("NTg" = "#F3D99E", "Q331K" = "#DBAEAF")) +
    
    # Error bars
    geom_errorbar(aes(ymin = .data[[group_col_name]] - SD,
                      ymax = .data[[group_col_name]] + SD),
                  width = 0.2) +
    
    # Facet
    facet_wrap(facet_formula, nrow=1, strip.position= "bottom", axes= "all") +
    
    # Graph titles
    labs(title = "Cytosol",
         x = "",
         y = "Relative expression (protein)",
         fill = x) + # Legend title
    
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0), labels = label_number(accuracy = 0.1))  # Setting both multiplier and add-on to 0
  
  # Define shapes for each DIFF value
  diff_shapes = c("F" = 21, "M" = 24)
  
  # Overlay individual data points (optional with different shapes for DIFF)
  p = p + geom_point(data=individual_data, aes(x = .data[[x]],
                                               y = .data[[individual_col_name]],
                                               shape = Sex),
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
  
  # Remove x-axis labels and ticks (OPTIONAL)
  p = p + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                strip.placement = "outside")  # Move facet labels below the plot
  
  # Add dashed line
  p = p + geom_hline(yintercept=1, linetype="dashed", color="black", size=0.3)
  
  # Control legend order
  p = p + guides(fill = guide_legend(order = 1),
                 shape = guide_legend(order = 2))
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_data(group_means, "Group_Mean", long_data, "norm_signal", "Genotype", "Protein")

plot

# Save plot
ggsave("/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/mouse_cytosol/wb_plot_all.png", plot=plot, width=12, height=3.5, dpi=300, bg="white")

# Export the test results to CSV files
write.csv(test_result, paste0(parent_filepath, "test_result.csv"), row.names=FALSE)

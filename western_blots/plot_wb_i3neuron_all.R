
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/i3neuron_synapse/"
filename = "total_allproteins.xlsx"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Define Replicate as a factor with levels
empiria_data$Replicate = factor(empiria_data$Replicate, levels = c("WT", "Q331K"))

# Reshape data into long format
long_data = empiria_data %>%
  pivot_longer(-Replicate, names_to = "Protein", values_to = "Signal") %>%
  arrange(Protein)

long_data$Protein = factor(long_data$Protein, levels = c("Synaptophysin", "Syntaxin-1A", "VAMP2", "SNAP-25", "UNC13A", "STXBP1", "PSD-95", "Homer-1"))

# Find the maximum y-values -> this is for dynamic annotation bars in the plots
max_y_per_protein = aggregate(Signal ~ Protein, data=long_data, max)
# Rename column
max_y_per_protein = max_y_per_protein %>%
  rename(max_y = Signal)


# Find group means --------------------------------------------------------

# Calculate mean for each group
group_means = long_data %>%
  group_by(Replicate, Protein) %>%
  summarise(
    N = n(), # Count the number of samples in each group
    Group_Mean = mean(`Signal`), # Calculate mean for each group
    .groups = "drop"
    )

# Find sd using normalised (to control) values
sd_values = long_data %>%
  group_by(Replicate, Protein) %>%
  summarise(
    SD = sd(Signal), .groups = "drop"
    )

# Merge `sd_values` with `group_means`
group_means = group_means %>%
  left_join(sd_values, by = c("Replicate", "Protein"))


# Statistics --------------------------------------------------------------

# Create function to test normality, equal variance, and perform appropriate statistical test
perform_test = function(data) {
  
  # Extract values grouped by Replicate using Signal
  # Test for normality using Shapiro-Wilk test
  normality_result = by(data$Signal, data$Replicate, shapiro.test)
  
  wt_normal = normality_result$WT$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                 # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
  q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                       # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
  
  if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
    # If data is not normally distributed, perform Mann-Whitney test
    message = "Data not normally distributed; perform Mann-Whitney test"
    test_result = wilcox.test(Signal ~ Replicate, data = data) # `Signal ~ Replicate` specifies that you want to compare values in the `Signal` column grouped by `Replicate`
  } else {
    # Test for equal variances using F-test
    variance_test = var.test(Signal ~ Replicate, data = data)
    var_equal = variance_test$p.value > 0.05 # If p-value is greater than 0.05, it assumes both replicates have equal variances
    
    if (!var_equal) { # If `var_equal` variable is set to `FALSE`
      # If variances are not equal, perform t-test with var.equal = FALSE
      message = "Data does not have equal variance; perform Welch's t-test"
      test_result = t.test(Signal ~ Replicate, data = data, var.equal = FALSE)
    } else {
      # If variances are equal, perform t-test with var.equal = TRUE
      message = "Data has equal variance; perform Student's t-test"
      test_result = t.test(Signal ~ Replicate, data = data, var.equal = TRUE)
    }
  }
  
  # Return both test results and the message
  return(list(test_result = test_result, message = message))
}

# Perform tests for each protein
# group_by and summarise are used to iterate over groups of data
test_results <- long_data %>%
  # Group by Protein column, which createas a subset of the dataset, where each subset corresponds to a unique protein
  group_by(Protein) %>%
  # Apply the perform_test function to subset of data
  # cur_data() refers to the current group's data subset
  summarise(
    p_value = perform_test(cur_data())$test_result$p.value,
    method = perform_test(cur_data())$test_result$method,
    alternative = perform_test(cur_data())$test_result$alternative,
    .groups = "drop"
  )

# Add significance stars
test_results <- test_results %>%
  mutate(
    Stars = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
  
# Merge the max y-values per protein with the test_results data frame
test_results = merge(test_results, max_y_per_protein, by = "Protein")

# Replace max_y with NA if not significant
test_results$max_y = ifelse(test_results$Stars == "", NA, test_results$max_y)

  
# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for facet bar plots
my_theme_facet = function() {
  theme_minimal() +
    theme(
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          panel.spacing = unit(0.25, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 10),  # Increase facet title size and make it bold
          axis.title.y = element_text(margin = margin(r = 10), # Adjust y-axis title position
                                      size = 10), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Calculate the maximum y value to set upper axis limit
max_y_value = max(group_means$Group_Mean + group_means$SD)
upper_limit = max_y_value * 1.25  # 25% buffer above the max value

# Define custom labeller function to add "Total " to facet titles
my_labeller = as_labeller(function(x) paste("Total", "\n", x))

plot = ggplot(group_means, aes_string(
                                x = "Replicate",
                                y = "Group_Mean",
                                fill = "Replicate"
                                )) +
  # Bar plot
  geom_col(position = position_dodge(0.5),
          width = 0.6,
          color = "black") +
  scale_fill_manual(values = c("WT" = "grey40", "Q331K" = "grey88")) +
  # Error bars
  geom_errorbar(aes(ymin=Group_Mean - SD, ymax=Group_Mean + SD),
                width=0.2,
                position=position_dodge(0.5)) +
  # Facet
  facet_wrap(~Protein, nrow=1, strip.position = "bottom") +
  # Graph titles
  labs(x="",
       y="Relative expression (total protein)",
       fill="Genotype") +
  # Plot appearance
  my_theme_facet() +
  scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0)) +  # Setting both multiplier and add-on to 0
  # Overlay individual data points
  geom_point(data = long_data, aes_string(x = "Replicate",
                                          y = "Signal"),
             position = position_dodge(0.5), size = 1.5) +
  # Significance stars
  geom_text(data = test_results, aes(label = Stars, x = 1.5, y = max_y + 0.2),
            position = position_dodge(width = 0.5), inherit.aes = FALSE, vjust = -0.5,
            size = 7) +  # Adjust size here
  # Significance lines
  geom_segment(data = test_results, aes(x = 1, xend = 2, y = max_y + 0.4, yend = max_y + 0.4),
                       linetype = "solid", color = "black", position = position_dodge(width = 0.5), inherit.aes = FALSE) +
  # Remove x-axis labels and ticks
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside") +  # Move facet labels below the plot
  # Add dashed line
  geom_hline(yintercept=1, linetype="dashed", color="black", size=0.3)

print(plot)

# Export plot
# Open a PNG file to save the plot
#
# For plot title with 1 line:
# width=825, height=1335
#
# For plot title with 2 lines:
# width=825, height=1390
png(paste0(parent_filepath, "wb_plot_all.png"), width=3500, height=1000, res=300)

# Create a plot
plot

# Close the device
dev.off()

# Export test results
# Function to extract p-value, method, alternative hypothesis, and sample sizes per group from test results
extract_test_results = function(test_results, data) {
  results_df = data.frame(
    Protein = character(),
    p_value = numeric(),
    method = character(),
    alternative = character(),
    WT_sample_size = integer(),  # Sample size for WT genotype
    Q331K_sample_size = integer(),  # Sample size for Q331K genotype
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(test_results))) {
    # Select current protein
    protein = test_results$Protein[i]
    
    # Subset data for current protein
    protein_data = subset(data, Protein == protein)
    
    # Extract the sample size for each genotype
    wt_sample_count <- nrow(subset(protein_data, Replicate == "WT"))
    q331k_sample_count <- nrow(subset(protein_data, Replicate == "Q331K"))
    
    # Add results for the current protein
    new_row <- data.frame(
      Protein = protein,
      p_value = round(test_results$p_value[i], 4),  # Round p-value to 4 decimal places
      method = test_results$method[i],  # Statistical test method
      alternative = test_results$alternative[i],  # Alternative hypothesis
      WT_sample_size = wt_sample_count,  # Sample size for WT genotype
      Q331K_sample_size = q331k_sample_count  # Sample size for Q331K genotype
    )
    # Combine
    results_df = rbind(results_df, new_row)
  }
  return(results_df)
}

# Extract results
results_df = extract_test_results(test_results, long_data)

# Define file path for saving CSVs
csv_path = paste0(parent_filepath, "test_result.csv")

# Export the test results to CSV files
write.csv(results_df, csv_path, row.names=FALSE)

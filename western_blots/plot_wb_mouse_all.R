
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/mouse_synaptosome/"
filename = "synaptosome_allproteins.xlsx"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Define Replicate as a factor with levels
empiria_data$Replicate = factor(empiria_data$Replicate, levels = c("NTg", "Q331K"))

# Reshape data into long format
long_data = empiria_data %>%
  pivot_longer(-c(Replicate, Sex), names_to = "Protein", values_to = "Signal") %>%
  arrange(Protein)

long_data$Protein = factor(long_data$Protein, levels = c("Synaptophysin", "Syntaxin-1A", "VAMP2", "SNAP-25", "Munc13-1", "Munc18-1", "PSD-95", "Homer-1"))

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
  
  NTg_normal = normality_result$NTg$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                 # `NTg_normal` is set to `TRUE` if p-value is greater than 0.05
  q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
                                                       # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
  
  if (!NTg_normal || !q331k_normal) { # If either variable is set to `FALSE`
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

# Filter for only complete data
long_data_filtered = long_data %>%
  filter(Signal != 0)

# Perform tests for each protein
test_results <- long_data_filtered %>%
  group_by(Protein) %>%
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
      plot.title = element_text(face = "bold", hjust = 0.5),# Adjust plot title
      axis.title.y = element_text(margin = margin(r = 15), # Adjust y-axis title position
                                  size = 12), # Adjust y-axis title size
      axis.text.x = element_text(size = 10), # Increase x-axis text size
      axis.text.y = element_text(size = 10), # Increase y-axis text size
      # Facet-specific
      panel.spacing = unit(0.5, "lines"),  # Adjust space between facet panels
      strip.text = element_text(size = 10)  # Increase facet title size
    ) 
}

# Calculate the maximum y value to set upper axis limit
max_y_value = max(group_means$Group_Mean + group_means$SD)
upper_limit = max_y_value * 1.25  # 25% buffer above the max value

# Define custom labeller function to add "Total " to facet titles
my_labeller = as_labeller(function(x) paste("Synaptosomal", "\n", x))

plot = ggplot(group_means, aes_string(
                                x = "Replicate",
                                y = "Group_Mean",
                                fill = "Replicate"
                                )) +
  # Bar plot
  geom_col(position = position_dodge(0.5),
          width = 0.8,
          color = "black") +
  scale_fill_manual(values = c("NTg" = "#F3D99E", "Q331K" = "#DBAEAF")) +
  # Error bars
  geom_errorbar(aes(ymin=Group_Mean - SD, ymax=Group_Mean + SD),
                width=0.2,
                position=position_dodge(0.5)) +
  # Facet
  facet_wrap(~Protein, nrow=1, strip.position = "bottom", axes="all") +
  # Graph titles
  labs(x="",
       y="Relative expression (protein)",
       fill="Genotype",
       title="Synaptosome") +
  # Plot appearance
  my_theme_facet() +
  scale_y_continuous(limits = c(0, 1.75), expand = c(0, 0)) +  # Setting both multiplier and add-on to 0
  # Overlay individual data points
  geom_point(data = long_data, aes_string(x = "Replicate",
                                          y = "Signal",
                                          shape = "Sex",
                                          group = "Replicate"),
             position = position_dodge(0.5), size = 1.9) +
  # Significance stars
  geom_text(data = test_results, aes(label = Stars, x = 1.5, y = (max_y + 0.065*upper_limit)),
            position = position_dodge(width = 0.5), inherit.aes = FALSE, vjust = -0.5,
            size = 6) +  # Adjust size here
  # Significance lines
  geom_segment(data = test_results, aes(x = 1, xend = 2, y = (max_y + 0.1*upper_limit), yend = (max_y + 0.1*upper_limit)),
               linetype = "solid", color = "black", position = position_dodge(width = 0.5), inherit.aes = FALSE) +
  # Remove x-axis labels
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside") +  # Move facet labels below the plot
  # Add dashed line
  geom_hline(yintercept=1, linetype="dashed", color="black", size=0.3) +
  # Control legend order
  guides(fill = guide_legend(order = 1, override.aes = list(shape = NA)),  # Genotype legend first, removing circle from fill legend
         shape = guide_legend(order = 2))  # Sex legend second

print(plot)

# Save plot
ggsave("/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/mouse_synaptosome/wb_plot_all.png", plot=plot, width=12, height=3.5, dpi=300, bg="white")

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
    ntg_sample_count <- nrow(subset(protein_data, Replicate == "NTg"))
    q331k_sample_count <- nrow(subset(protein_data, Replicate == "Q331K"))
    
    # Add results for the current protein
    new_row <- data.frame(
      Protein = protein,
      p_value = round(test_results$p_value[i], 4),  # Round p-value to 4 decimal places
      method = test_results$method[i],  # Statistical test method
      alternative = test_results$alternative[i],  # Alternative hypothesis
      NTg_sample_size = ntg_sample_count,  # Sample size for WT genotype
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


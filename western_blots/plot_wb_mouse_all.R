
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/mouse_synaptosome/"
filename = "synaptosome_allproteins.xlsx"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Define Replicate as a factor with levels
empiria_data$Replicate = factor(empiria_data$Replicate, levels = c("NTg", "Q331K"))

# Reshape data into long format
long_data = empiria_data %>%
  pivot_longer(-Replicate, names_to = "Protein", values_to = "Signal") %>%
  arrange(Protein)

long_data$Protein = factor(long_data$Protein, levels = c("Synaptophysin", "Syntaxin 1A", "VAMP2", "SNAP25", "Munc13-1", "Munc18-1", "PSD-95", "Homer"))


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

# Perform tests for each protein
test_results <- long_data %>%
  group_by(Protein) %>%
  summarise(
    p_value = perform_test(cur_data())$test_result$p.value,
    message = perform_test(cur_data())$message,
    .groups = "drop"
  )

# Add significance stars
test_results <- test_results %>%
  mutate(
    Stars = case_when(
      p_value <= 0.0001 ~ "****",
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Find the maximum y-values -> this is for dynamic annotation bars in the plots
test_results = test_results %>%
  mutate(
    max_y = if_else(Stars != "", max(group_means$Group_Mean + group_means$SD), NA_real_)
  )


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
upper_limit = max_y_value * 1.5  # 50% buffer above the max value

# Define custom labeller function to add "Total " to facet titles
my_labeller = as_labeller(function(x) paste("Synaptosomal", "\n", x))

plot = ggplot(group_means, aes_string(
                                x = "Replicate",
                                y = "Group_Mean",
                                fill = "Replicate"
                                )) +
  # Bar plot
  geom_col(position = position_dodge(0.5),
          width = 0.6,
          color = "black") +
  scale_fill_manual(values = c("NTg" = "grey40", "Q331K" = "grey88")) +
  # Error bars
  geom_errorbar(aes(ymin=Group_Mean - SD, ymax=Group_Mean + SD),
                width=0.2,
                position=position_dodge(0.5)) +
  # Facet
  facet_wrap(~Protein, nrow=1, strip.position = "bottom") +
  # Graph titles
  labs(x="",
       y="Relative expression (synaptosomal protein)",
       fill="Genotype") +
  # Plot appearance
  my_theme_facet() +
  scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0)) +  # Setting both multiplier and add-on to 0
  # Overlay individual data points
  geom_point(data = long_data, aes_string(x = "Replicate",
                                            y = "Signal"),
             position = position_dodge(0.5), size = 1.5) +
  # Significance stars
  geom_text(data = test_results, aes(label = Stars, x = 1.5, y = max_y * 1),
            position = position_dodge(width = 0.5), inherit.aes = FALSE, vjust = -0.5,
            size = 7) +  # Adjust size here
  # Significance lines
  geom_segment(data = test_results, aes(x = 1, xend = 2, y = max_y * 1.05, yend = max_y * 1.05),
                       linetype = "solid", color = "black", position = position_dodge(width = 0.5), inherit.aes = FALSE) +
  # Remove x-axis labels
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

# Print message
print(message)

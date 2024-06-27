
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/i3neuron_synapse/synper_extraction_practice/"
relative_filepath = "unc/"
filename = "empiria_unc.xlsx"
protein_name = "UNC13A"

colour1 = "grey40"
colour2 = "grey88"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
empiria_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Prepare data frame
empiria_data = empiria_data[!is.na(empiria_data[, 4]), ] # Directly removes any rows with any `NA` values in the fourth column (Normalised Signal)
                                                         # Only keep the rows where fourth column is not NA
colnames(empiria_data) = as.character(empiria_data[1,]) # Make the first row the column names
empiria_data = empiria_data[-1, ] # Remove first row

# Split the Name column
empiria_data = empiria_data %>%
  separate(Name, into = c("Replicate", "Fraction"), sep = " ")

# Define `Normalized Signal` variable as numeric
empiria_data$`Normalized Signal` = as.numeric(empiria_data$`Normalized Signal`)

# Define Replicate and Fraction as a factor with levels
empiria_data$Replicate = factor(empiria_data$Replicate, levels = c("WT", "Q331K"))
empiria_data$Fraction = factor(empiria_data$Fraction, levels = c("Lysate", "Cytosol", "Synaptosome"))


# # Find group means --------------------------------------------------------
# 
# # Calculate mean for each group
# group_means = empiria_data %>%
#   group_by(Replicate) %>%
#   summarise(
#     N = n(), # Count the number of samples in each group
#     Group_Mean = mean(`Normalized Signal`), # Calculate mean for each group
#     )
#  
# 
# # Normalise samples to control mean ---------------------------------------
# 
# # Find the mean of the control group
# control_mean = group_means %>%
#   filter(Replicate == "WT") %>% # Where `Replicate` equals `WT`
#   pull(Group_Mean) # Extracts `Group_Mean` column (in this case, a single value)
# 
# # Divide each group mean by the control mean
# group_means = group_means %>%
#   mutate(Normalised_Group_Mean = Group_Mean / control_mean)
# 
# # Divide each sample by the control mean
# normalised_to_control = empiria_data %>%
#   select(Replicate, `Normalized Signal`) %>% # Select the `Replicate` and `Normalized Signal` columns
#   mutate(Normalised_to_Control = `Normalized Signal` / control_mean) %>%
#   arrange(Replicate)
# 
# # Export .csv
# write.csv(normalised_to_control, paste0(parent_filepath, relative_filepath, protein_name, "_normalised_to_control.csv"), row.names=FALSE)
# 
# # Find sd using normalised (to control) values
# sd_values = normalised_to_control %>%
#   group_by(Replicate) %>%
#   summarise(
#     SD = sd(Normalised_to_Control)
#     )
# 
# # Merge `sd_values` with `group_means`
# group_means = group_means %>%
#   inner_join(sd_values, by = "Replicate")
# 
# 
# # Statistics --------------------------------------------------------------
# 
# # Create function to test normality, equal variance, and perform appropriate statistical test
# perform_test = function(data) {
#   
#   # Extract values grouped by Replicate using Normalised_to_Control
#   # Test for normality using Shapiro-Wilk test
#   normality_result = by(data$Normalised_to_Control, data$Replicate, shapiro.test)
#   
#   wt_normal = normality_result$WT$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
#                                                  # `wt_normal` is set to `TRUE` if p-value is greater than 0.05
#   q331k_normal = normality_result$Q331K$p.value > 0.05 # If p-value is greater than 0.05, the data is considered normally distributed
#                                                        # `qk_normal` is set to `TRUE` if p-value is greater than 0.05
#   
#   if (!wt_normal || !q331k_normal) { # If either variable is set to `FALSE`
#     # If data is not normally distributed, perform Mann-Whitney test
#     message = "Data not normally distributed; perform Mann-Whitney test"
#     test_result = wilcox.test(Normalised_to_Control ~ Replicate, data = data) # `Normalised_to_Control ~ Replicate` specifies that you want to compare values in the `Normalised_to_Control` column grouped by `Replicate`
#   } else {
#     # Test for equal variances using F-test
#     variance_test = var.test(Normalised_to_Control ~ Replicate, data = data)
#     var_equal = variance_test$p.value > 0.05 # If p-value is greater than 0.05, it assumes both replicates have equal variances
#     
#     if (!var_equal) { # If `var_equal` variable is set to `FALSE`
#       # If variances are not equal, perform t-test with var.equal = FALSE
#       message = "Data does not have equal variance; perform Welch's t-test"
#       test_result = t.test(Normalised_to_Control ~ Replicate, data = data, var.equal = FALSE)
#     } else {
#       # If variances are equal, perform t-test with var.equal = TRUE
#       message = "Data has equal variance; perform Student's t-test"
#       test_result = t.test(Normalised_to_Control ~ Replicate, data = data, var.equal = TRUE)
#     }
#   }
#   
#   # Return both test results and the message
#   return(list(test_result = test_result, message = message))
# }
# 
# # Perform t-test
# result = perform_test(normalised_to_control)
# 
# # Accessing the test results
# test_result = result$test_result
# 
# # Accessing the message
# message = result$message
# 
# # Create function that converts p-values to signficance stars and generates a dataframe of annotations
# prepare_annotations = function(test_result) {
#   # Convert p-values to stars based on traditional significance levels
#   convert_p_to_stars = function(p_value) {
#     if (p_value <= 0.0001) {
#       return("****")
#     } else if (p_value <= 0.001) {
#       return("***")
#     } else if (p_value <= 0.01) {
#       return("**")
#     } else if (p_value <= 0.05) {
#       return("*")
#     } else {
#       return("")  # Not significant
#     }
#   }
#   
#   # Extract p-value from test results
#   p_value = test_result$p.value
#   
#   # Convert p-value to stars
#   stars = convert_p_to_stars(p_value)
#   
#   # Create a dataframe `annotations`
#   annotations = data.frame(
#     p_value = p_value,
#     Stars = stars,
#     stringsAsFactors = FALSE
#   )
#   
#   # Filter out non-significant annotations
#   annotations = annotations[annotations$Stars != "", ]
#   
#   return(annotations)
# }
# 
# # Create annotations 
# annotation = prepare_annotations(test_result)
# 
# # Find the maximum y-values -> this is for dynamic annotation bars in the plots
# annotation = annotation %>%
#   mutate(
#     max_y = max(group_means$Normalised_Group_Mean + group_means$SD)
#   )


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for facet bar plots
my_theme_facet = function() {
  theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          panel.spacing = unit(1, "lines"),  # Adjust space between facet panels
          strip.text = element_text(size = 12, face = "bold"),  # Increase facet title size and make it bold
          axis.title.y = element_text(margin = margin(r = 10), # Adjust y-axis title position
                                      size = 12), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Create a function that takes two dataframes and column names to generate multiple bar plots with overlayed data points
plot_normalised = function(empiria_data, x = "Fraction", y = "`Normalized Signal`") {
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(empiria_data$`Normalized Signal`, na.rm = TRUE)
  upper_limit = max_y_value * 1.1 # 10% buffer above the max value
  
  # Create the bar plot
  p = ggplot(empiria_data, aes_string(x = x,
                                    y = y,
                                    fill = "Replicate")) +
    # Bar plot
    geom_col(position = position_dodge(0.9),
             width = 0.6,
             color = "black") +
    scale_fill_manual(values = c("WT" = colour1, "Q331K" = colour2)) +
    
    # Facet by Replicate
    facet_wrap(~Replicate) +
    
    # Graph titles
    labs(title = protein_name,
         x = "",
         y = "Normalised Signal",
         fill = x) +

    # Plot appearance
    my_theme_facet() +
    scale_y_continuous(limits = c(0, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
  
  # Print the plot
  return(p)
}

# Make plot
plot = plot_normalised(empiria_data)

plot

# Export plot
# Open a PNG file to save the plot
png(paste0(parent_filepath, relative_filepath, protein_name, "_enrichment.png"), width=1735, height=1335, res=300)

# Create a plot
plot

# Close the device
dev.off()

# # Print message
# print(message)

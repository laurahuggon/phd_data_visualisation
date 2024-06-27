
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)


# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/wb/i3neuron_synapse/"
relative_filepath = "synper_extraction_practice/"
filename = "synper_extraction_practice.xlsx"


# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
protein_conc_data = read_excel(full_filename)


# Data processing ---------------------------------------------------------

# Calculate the mean absorbance for standards and samples
protein_conc_data = protein_conc_data %>%
  # Calculate mean of all columns starting with "Absorbance", ignoring `NA` values
  mutate(Avg_Absorbance = rowMeans(select(., starts_with("Absorbance")), na.rm = TRUE))

# Filter out the standards data
standards = protein_conc_data %>%
  filter(!is.na(Concentration))


# Data visualisation 1 ----------------------------------------------------

# Create custom ggplot2 theme
my_theme = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold",
                                    hjust = 0.5), # Adjust plot title
          axis.title.y = element_text(margin = margin(r = 14), # Adjust y-axis title position
                                      size = 13), # Adjust y-axis title size
          axis.title.x = element_text(margin = margin(t = 14), # Adjust x-axis title position
                                      size = 13), # Adjust x-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Plot the standard curve using the `Concentration` and `Avg_Absorbance` values for the standards.
plot = ggplot(data = protein_conc_data, aes(x = Concentration, y = Avg_Absorbance)) +
  # Linear regression line through the standards
  geom_smooth(data = filter(protein_conc_data, Data_Type == "Standard"),
              method = "lm",
              color = "black",
              size = 0.5) +
  
  # Scatter plot for standards
  geom_point(data = filter(protein_conc_data, Data_Type == "Standard"),
             color = "black") +
  
  # Graph titles
  labs(title = "Protein standard curve",
       x = "Protein concentration (mg/mL)",
       y = "Absorbance (750 nm)") +
  
  # Plot appearance
  my_theme() 


# Linear regression model -------------------------------------------------

# Fit linear model
model = lm(Avg_Absorbance ~ Concentration, data = standards) # Absorbance is explained by concentration

# Extracting model coefficients that describe the equation of the line
intercept = coef(model)[1]
slope = coef(model)[2]

# Create equation string
# `sprintf()` is used to format the equation string, incorporating the coefficients into the linear formula y = mx + b
# where m is the slope and b is the intercept
equation = sprintf("y = %.3fx + %.3f", slope, intercept)

# Interpolate concentrations for samples ----------------------------------

# Filter out sample data
samples = protein_conc_data %>%
  filter(is.na(Concentration))

# Interpolate concentrations, rounded to 2 dp
samples$Interpolated_Concentration = round(((samples$Avg_Absorbance - intercept) / slope), 2) # Estimate concentration using the equation derived from the regression line

# Update original dataframe
# If `Data_Type` is "Sample", replace `Concentration` with values from `samples$Interpolated_Concentration`
# The `match()` function aligns indices of the rows labelled as "Sample" in `protein_conc_data` with the rows in `samples` dataframe
# `row_number()` returns a vector of the row indicies in `protein_conc_data` and `which(Data_Type == "Sample"`) returns the indices of just the "Sample" rows
protein_conc_data = protein_conc_data %>%
  mutate(Concentration = if_else(Data_Type == "Sample",
                                 samples$Interpolated_Concentration[match(row_number(), which(Data_Type == "Sample"))], 
                                 Concentration))


# Data visualisation 2 ----------------------------------------------------

# Select the `Avg_Absorbance` and `Interpolated_Concentration` columns from `samples` data to be used in the plot
samples = samples %>%
  select(Avg_Absorbance, Interpolated_Concentration)

# Calculate positions
samples$y_position = 0.24 - seq(0, by = 0.02, length.out = nrow(samples))  # Adjust spacing as needed

# Modify plot
plot = plot +
  # Scatter plot for samples
  geom_point(data = filter(protein_conc_data, Data_Type == "Sample"),
             color = "red") +
  
  # Equation text
  annotate("text", x = 0.15, y = 0.3, label = equation, size = 4, color = "black", hjust = 0, vjust = 0) +
  
  # Column titles
  annotate("text", x = 1.25, y = 0.26, label = "Abs", hjust = 0.5, vjust = 1, size = 4) +
  annotate("text", x = 1.425, y = 0.26, label = "Conc", hjust = 0.5, vjust = 1, size = 4) +
  
  # Annotations for absorbance and concentration values
  geom_text(data = samples, aes(x = 1.25, y = y_position, label = sprintf("%.3f", Avg_Absorbance)), hjust = 0.5, vjust = 1, size = 4, color = "black") +
  geom_text(data = samples, aes(x = 1.425, y = y_position, label = sprintf("%.2f", Interpolated_Concentration)), hjust = 0.5, vjust = 1, size = 4, color = "black")

plot

# Export plot
# Open a PNG file to save the plot
png(paste0(parent_filepath, relative_filepath, "dc_standard_curve.png"), width=2221.8, height=1390, res=300)

# Create a plot
plot

# Close the device
dev.off()

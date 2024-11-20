# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggrepel)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King\'sCollegeLondon/phd/lab/wb/i3neuron_synapse/"
relative_filepath = "bca_test/"
filename = "bca_abs_noconc.xlsx"

# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
protein_conc_data = read_excel(full_filename)


# Data processing ---------------------------------------------------------

# Select columns/rows needed
protein_conc_data = protein_conc_data %>%
  select(Data_Type, Concentration, Absorbance_1, Absorbance_2)
protein_conc_data = protein_conc_data[-1, ]

# Ensure concentration is numeric
protein_conc_data = protein_conc_data %>%
  mutate(Concentration = as.numeric(Concentration))

# Calculate the mean absorbance for standards and samples
protein_conc_data = protein_conc_data %>%
  # Ensure columns are numeric
  mutate(across(starts_with("Absorbance"), as.numeric)) %>%
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
  # Concentration cutoff line
  geom_rect(
    aes(xmin=-Inf, xmax=125, ymin=-Inf, ymax=Inf),
    fill="grey", alpha=0.05
  ) +
  geom_vline(xintercept=125, linetype="dashed", color="black", size=0.3) +
  
  # Quadratic regression line
  geom_smooth(data = filter(protein_conc_data, Data_Type == "Standard"), # Only use standard data points for the regression
              formula = y ~ poly(x, 2), # Second-degree polynomial (quadratic regression)
              method = "lm", # A linear model is used to fit the polynomial
              color = "black",
              size = 0.5) +
  
  # Scatter plot for standards
  geom_point(data = filter(protein_conc_data, Data_Type == "Standard"),
             color = "black", shape=15) +
  
  # Graph titles
  labs(title = "Protein standard curve",
       x = "Protein concentration (μg/mL)",
       y = "Absorbance (562 nm)") +
  
  # Plot appearance
  my_theme() 

plot

# Quadratic regression model -------------------------------------------------

# Fit quadratic model using the standards, modelling Avg_Absorbance as a quadratic function of Concentration
model = lm(Avg_Absorbance ~ poly(Concentration, 2, raw = TRUE), data = standards) # Second-degree polynomial with raw (unorthogonalise) coefficients

# Extracting model coefficients for the quadratic equation (y = ax^2 + bx + c)
intercept = coef(model)[1] # c
coef_x = coef(model)[2] # b
coef_x2 = coef(model)[3] # a

# Create equation string for a quadratic formula
equation = sprintf("y = %.3fx^2 + %.3fx + %.3f", coef_x2, coef_x, intercept)

# Interpolate concentrations for samples ----------------------------------

# Filter out sample data
samples = protein_conc_data %>%
  filter(is.na(Concentration))

# Define a function to solve the quadratic equation for concentration
solve_concentration <- function(absorbance, intercept, coef_x, coef_x2) {
  # Solves ax^2 + bx + c = 0
  a = coef_x2
  b = coef_x
  c = intercept - absorbance
  discriminant = b^2 - 4 * a * c
  if (discriminant >= 0) {
    # Only return the positive root - only the positive root is relevant for concentration values
    conc = (-b + sqrt(discriminant)) / (2 * a)
    return(conc)
  } else {
    return(NA)
  }
}

# Interpolate concentrations for samples and round to 0 decimal places
samples$Interpolated_Concentration = round(
  # Loop through each sample's Avg_Absorbance and apply solve_concentration
  sapply(samples$Avg_Absorbance, solve_concentration, intercept, coef_x, coef_x2), 
  digits = 0
)

# Update original dataframe for samples only
protein_conc_data = protein_conc_data %>%
  mutate(Concentration = if_else(Data_Type == "Sample",
                                 # Find the corresponding sample row (in original df) for which interpolated concentrations values should be assigned
                                 # row_number() generates a vector of row numbers for the entire dataframe (protein_conc_data)
                                 # which() returns the indices (row numbers) of all the "Sample" rows
                                 # match() finds whcih rows in protein_conc_data correspond to sample rows by comparing the row numbers
                                  # from rownumber() with the row numbers obtained from which()
                                 samples$Interpolated_Concentration[match(row_number(), which(Data_Type == "Sample"))],
                                 # For rows that are not samples, the original Concentration is kept
                                 Concentration))


# Data visualisation 2 ----------------------------------------------------

# Select the `Avg_Absorbance` and `Interpolated_Concentration` columns from `samples` data to be used in the plot
samples = samples %>%
  select(Avg_Absorbance, Interpolated_Concentration)

# Modify plot
plot = plot +
  # Scatter plot for samples
  geom_point(data = filter(protein_conc_data, Data_Type == "Sample"),
             color = "red") +
  
  # Equation text
  annotate("text", x = 1200, y = 0.2, label = equation, size = 4, color = "black", hjust = 0, vjust = 0) +
  
  # Absorbance values above each sample point
  geom_text_repel(data = filter(protein_conc_data, Data_Type == "Sample"),
                  aes(x = Concentration, y = Avg_Absorbance, label = sprintf("%.3f", Avg_Absorbance)),
                  nudge_y = 0.3, size = 3, color = "red", direction = "y", vjust = -0.5, segment.alpha = 0.25) +
  
  # Concentration values below each sample point (only if conc is >= 0.2)
  geom_text_repel(data = filter(protein_conc_data, Data_Type == "Sample", Concentration >= 125),
                  aes(x = Concentration, y = Avg_Absorbance, label = sprintf("%.0f", Concentration)),
                  nudge_y = -0.3, size = 3, color = "red", direction = "y", vjust = 4, segment.alpha = 0.25)

plot

# Export plot
# Open a PNG file to save the plot
png(paste0(parent_filepath, relative_filepath, "bca_standard_curve.png"), width=2221.8, height=1390, res=300)

# Create a plot
plot

# Close the device
dev.off()
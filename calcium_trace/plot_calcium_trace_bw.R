
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

# Define variables --------------------------------------------------------

parent_filepath = "/Users/laurahuggon/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y2/calcium_imaging/"
relative_filepath = "23.12.07/traces/"
filename = "LAURA_TDP43_iCortical_norm_neurobasal.xlsx"
genotype = "Q331K"


# Load data ---------------------------------------------------------------

full_filename = paste0(parent_filepath, relative_filepath, filename)
trace_data = read_excel(full_filename)


# Prepare data ------------------------------------------------------------

# Prepare data frame
trace_data = na.omit(trace_data) # Directly removes any rows with any `NA` values
colnames(trace_data) = as.character(trace_data[1,]) # Make the first row the column names
trace_data = trace_data[-1, ] # Remove first row
trace_data = trace_data[,-1] # Remove first column

# Convert all character columns to numeric
trace_data = trace_data %>% 
  mutate(across(where(is.character), as.numeric))

# Dynamically rename the columns from the 3rd column onwards
# Generates names starting with "Cell" followed by a sequence number
# The sequence starts from 1 up to the number of columns minus the first one (as this should not to be renamed)
# Creates a complete vector of new column names, including the unchanged one
new_column_names = c("Time [s]", paste("Cell", 1:(ncol(trace_data)-1)))

# Apply the new names to the dataframe
trace_data = setNames(trace_data, new_column_names)


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for bar plots
my_theme = function() {
  theme_minimal() +
    theme(axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          plot.title = element_text(face = "bold",
                                    hjust = 0.5), # Adjust plot title
          axis.title.y = element_text(margin = margin(r = 15), # Adjust y-axis title position
                                      size = 13), # Adjust y-axis title size
          axis.text.x = element_text(size = 10), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}


# Create a function that 
plot_trace = function(data) {
  # Make sure 'data' is a tibble, if not, convert it to one
  data = as_tibble(data)
  
  # Gather the cell columns into a long format suitable for plotting
  data_long = data %>%
    pivot_longer(
      cols = -`Time [s]`, # Exclude the time column from gathering
      names_to = "Cell", # New column for cell names
      values_to = "Intensity" # New column for intensity values
    )
  
  # Calculate the maximum y value to set upper axis limit
  max_y_value = max(data_long[["Intensity"]], na.rm = TRUE)
  upper_limit = max_y_value * 1.25  # 2.5% buffer above the max value
  
  # Calculate the minimum y value to set upper axis limit
  min_y_value = min(data_long[["Intensity"]], na.rm = TRUE)
  lower_limit = min_y_value * 0.975  # 2.5% buffer below the min value
  
  # Create line plot
  ggplot(data_long, aes(x = `Time [s]`, y = Intensity, group = Cell)) +
    geom_line(color = "black") +  # Add line geometries
   
    # Graph titles
     labs(title = genotype,
         x = "Time (s)",
         y = "Mean fluo-4 intensity") +
     
    # Plot appearance
    my_theme() +
    scale_y_continuous(limits = c(lower_limit, upper_limit), expand = c(0, 0))  # Setting both multiplier and add-on to 0
}

# Make plot
plot = plot_trace(trace_data)

plot

# Export plot
# Open a PNG file to save the plot
png(paste0(parent_filepath, relative_filepath, genotype, "_calcium_trace_bw.png"), width=2150, height=1335, res=300)

# Create a plot
plot

# Close the device
dev.off()


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(ggbeeswarm)


# Define variables --------------------------------------------------------

protein1_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/eq_diameter/COLOCPRE_EqDia_syp.csv"
protein2_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/eq_diameter/COLOCPRE_EqDia_stx.csv"
protein3_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/eq_diameter/COLOCPOST_EqDia_psd.csv"
protein4_filepath = "/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/eq_diameter/COLOCPOST_EqDia_hmr.csv"

protein1 = "Synaptophysin"
protein2 = "Syntaxin-1A"
protein3 = "PSD-95"
protein4 = "Homer-1"


# Load data ---------------------------------------------------------------

protein1_eq_diameter = read_csv(protein1_filepath)
protein2_eq_diameter = read_csv(protein2_filepath)
protein3_eq_diameter = read_csv(protein3_filepath)
protein4_eq_diameter = read_csv(protein4_filepath)


# Prepare data ------------------------------------------------------------

# Add protein name to each dataframe
add_protein = function(data, protein_name){
  data = data %>%
    mutate(
      Protein = protein_name
    )
  
  return(data)
}

protein1_eq_diameter = add_protein(protein1_eq_diameter, protein1)
protein2_eq_diameter = add_protein(protein2_eq_diameter, protein2)
protein3_eq_diameter = add_protein(protein3_eq_diameter, protein3)
protein4_eq_diameter = add_protein(protein4_eq_diameter, protein4)

# Combine dataframes
eq_diameter_df = rbind(protein1_eq_diameter, protein2_eq_diameter, protein3_eq_diameter, protein4_eq_diameter)

# Convert um to nm
eq_diameter_df = eq_diameter_df %>%
  mutate(EqDiameter = EqDiameter * 1000)

# Find stats for each protein
eq_diameter_stats = eq_diameter_df %>%
  group_by(Protein) %>%
  summarise(
    N = n(),
    EqDiameter_Mean = mean(EqDiameter),
    Max = max(EqDiameter),
    Min = min(EqDiameter)
  ) %>%
  rename(EqDiameter = EqDiameter_Mean)

# Set `Protein` as a factor with levels
eq_diameter_stats$Protein = factor(eq_diameter_stats$Protein, levels = c("Synaptophysin", "Syntaxin-1A", "PSD-95", "Homer-1"))


# Data visualisation ------------------------------------------------------

# Create custom ggplot2 theme for scatter plots 
my_theme = function() {
  theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),  # Add axis lines
          axis.ticks = element_line(colour = "black"),  # Add axis ticks
          #plot.title = element_text(face = "bold",
                                    #hjust = 0.5), # Adjust plot title
          axis.title.y = element_text(margin = margin(r = 21), # Adjust y-axis title position
                                      size = 13), # Adjust y-axis title size
          axis.text.x = element_text(size = 10, margin = margin(t = 8)), # Increase x-axis text size
          axis.text.y = element_text(size = 10) # Increase y-axis text size
    ) 
}

# Set colours for each protein
custom_colors = c(
  "Synaptophysin" = "seagreen4",
  "Syntaxin-1A" = "seagreen4",
  "PSD-95" = "firebrick2",
  "Homer-1" = "firebrick2"
)

# Create the initial plot with individual points (scatter)
plot = ggplot(eq_diameter_df, aes(x = Protein, y = EqDiameter, fill = Protein, color = Protein)) +
  geom_quasirandom(alpha=0.75, size=0.5, groupOnX=TRUE) +
  #geom_violin(width = 0.8, alpha = 0, color = "Black", fill = "white") +
  labs(x = "",
       y = "Equivalent Diameter (nm)") +
  
  # Mean line
  geom_crossbar(data = eq_diameter_stats, aes(x=Protein, y=EqDiameter, ymin=EqDiameter, ymax=EqDiameter),
                width=0.9, color="black", fatten=1) +
  
  # Annotate mean, max and min values
  geom_text(data = eq_diameter_stats, aes(x = Protein, y = EqDiameter, label = sprintf("%.0f nm", EqDiameter)), vjust=-0.8, color = "black", size = 3.5) +
  geom_text(data = eq_diameter_stats, aes(x = Protein, y = Max, label = sprintf("%.0f nm", Max)), vjust = -0.8, color = "black", size = 3.5) +
  geom_text(data = eq_diameter_stats, aes(x = Protein, y = Min, label = sprintf("%.0f nm", Min)), vjust = 1.8, color = "black", size = 3.5) +


# Plot appearance
  my_theme() +
  #scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors)

plot

# Export plot
# Open a PNG file to save the plot
png(paste0("/Users/k21224575/Library/CloudStorage/OneDrive-King'sCollegeLondon/phd/lab/imaging/isim/imaging_data_y1/syp_stx/analysis_nis_elements/eq_diameter/plot.png"), width=2000, height=1700, res=300)

# Create a plot
plot

# Close the device
dev.off()


# PhD Experiments: Data Processing & Data Visualisation

This repository contains all the R scripts used for data processing and visualisation for experiments performed during my PhD. The repository is organised into folders, each dedicated to a specific type of analysis or dataset. Below is a brief description of each folder and scripts.

## Repository Structure

### calcium_trace
Contains an R script for visualising calcium traces from fluorescent imaging data. The script takes normalised Fluo-4 intensity values of ROIs measured over time and plots the trace to show calcium fluctuations.

### neuron_quantification
Contains an R script for plotting the % of marker+ cells from fluorescent imaging data. Finds sample and group means, tests normality and variance to perform statistical testing, and plots the results on a bar plot.

### protein_standard_curve
Contains R scripts used to generate and analyse protein standard curves from protein concentration assays. The script takes raw absorbance values for standards and samples to plot a standard curve (either linear or quadratic regression) and extrapolates unknown protein concentrations.

### rna_quantification
Contains R scripts for plotting the quality control statistics for RNA sample characterisation. Generates 4 plots for each metric (concentration, A260/A280, A260/A230, RINe), faceted by fraction and *n* number. There are two scripts - one for the optimisation extractions and one for the experimental extractions.

### synaptic_protein_cellbody_intensity
Contains R scripts for plotting the intensity of synaptic protein localised to the cell body. The scripts take intensity values generated from iSIM super-resolution fluorescent images that have been analysed using NIS-Elements (Nikon) software. Finds sample (by *n* number and genotype) and group (by genotype) means, tests normality and variance to perform statistical testing, and plots the results on a bar plot.
* `plot_cellbody_intensity_allimages.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means).
* `plot_cellbody_intensity_colour.R`: individual data points (independent differentiations) on plots are coloured by staining-to-imaging times. Green indicates a short staining-to-imaging time and red indicates a long staining-to-imaging time relative to all samples in the experiment.
* `plot_cellbody_intensity_n1.R`: generates plots for experiments containing *n*=1. Does not perform statistical testing.
* `plot_cellbody_intensity_relativecolour.R`: individual data points (independent differentiations) on plots are coloured by relative staining-to-imaging times. Green indicates a short staining-to-imaging time and red indicates a long staining-to-imaging time for mutant samples relative to its paired wildtype sample.
* `plot_cellbody_intensity.R`: generates plots for experiments containing *n*=3. Performs statistical testing.

### synaptic_protein_global_intensity
Contains R scripts for plotting the intensity of synaptic protein localised across the whole cell. The scripts take intensity values generated from iSIM super-resolution fluorescent images that have been analysed using NIS-Elements (Nikon) software. Finds sample (by *n* number and genotype) and group (by genotype) means, tests normality and variance to perform statistical testing, and plots the results on a bar plot.
* `plot_global_intensity_allimages.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means).
* `plot_global_intensity_colour.R`: individual data points (independent differentiations) on plots are coloured by staining-to-imaging times. Green indicates a short staining-to-imaging time and red indicates a long staining-to-imaging time relative to all samples in the experiment.
* `plot_global_intensity_n1.R`: generates plots for experiments containing *n*=1. Does not perform statistical testing.
* `plot_global_intensity_relativecolour.R`: individual data points (independent differentiations) on plots are coloured by relative staining-to-imaging times. Green indicates a short staining-to-imaging time and red indicates a long staining-to-imaging time for mutant samples relative to its paired wildtype sample.
* `plot_global_intensity.R`: generates plots for experiments containing *n*=3. Performs statistical testing.
* **tests**
  * `plot_global_intensity_median_matching.R`: to account for differences in staining-to-imaging times affecting fluorescent intensity, this script attempts to match medians between unaffected and affected samples.
  * `plot_global_intensity_normalised_to_control.R`: to account for differences in staining-to-imaging times affecting fluorescent intensity, this script normalises the mutant sample means to its paired wildtype sample mean.

### synaptic_protein_puncta_intensity
Contains R scripts for plotting the intensity of synaptic protein localised to synaptic puncta. The scripts take intensity values generated from iSIM super-resolution fluorescent images that have been analysed using NIS-Elements (Nikon) software. Finds sample (by *n* number and genotype) and group (by genotype) means, tests normality and variance to perform statistical testing, and plots the results on a bar plot.
* `plot_puncta_intensity_allimages.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means).
* `plot_puncta_intensity_n1.R`: generates plots for experiments containing *n*=1. Does not perform statistical testing.
* `plot_puncta_intensity.R`: generates plots for experiments containing *n*=3. Performs statistical testing.

### synaptic_puncta_morphology
Contains R scripts for plotting puncta morphology metrics (% colocalised puncta, the density of puncta, and the volume of puncta). The scripts take intensity values generated from iSIM super-resolution fluorescent images that have been analysed using NIS-Elements (Nikon) software. Finds sample (by *n* number and genotype) and group (by genotype) means, tests normality and variance to perform statistical testing, and plots the results on a bar plot for each DIV.
* `plot_eq_diameter.R`: generates a plot of the equivalent diameter of all puncta detected in one image.
* `plot_morphology_allimages.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means).
* `plot_morphology_coloc_allimages_div28.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means). Only colocalised puncta have been analysed. Includes DIV28 data.
* `plot_morphology_coloc_allimages.R`: performs statistical testing and generates plots using all the images in the dataset (rather than using the sample means). Only colocalised puncta have been analysed.
* `plot_morphology_coloc_n1.R`: generates plots for experiments containing *n*=1. Does not perform statistical testing. Only colocalised puncta have been analysed.
* `plot_morphology_coloc.R`: generates plots for experiments containing *n*=3. Performs statistical testing. Only colocalised puncta have been analysed.
* `plot_morphology_n1.R`: generates plots for experiments containing *n*=1. Does not perform statistical testing. All detected have been analysed.
* `plot_morphology.R`: generates plots for experiments containing *n*=3. Performs statistical testing. All detected have been analysed.

### western_blots
Contains R scripts for plotting normalised band intensity values from western blots that have been analysed using Empiria Studio (LI-COR) software. Finds sample and group means, normalises the mutant group mean to the wildtype group mean, tests normality and variance to perform statistical testing, and plots the results on a bar plot.
* `plot_wb_i3neuron_nonormalisation`: generates plots for experiments performed on i3Neuron whole-cell lysate samples. Normalised signal is plotted on the y-axis instead of fold change.
* `plot_wb_i3neuron_synaptosome_validation`: generates plots for validation experiments containing *n*=1. Normalised signal is plotted on the y-axis instead of fold change.
* `plot_wb_i3neuron`: generates plots for experiments performed on i3Neuron whole-cell lysate samples.
* `plot_wb_mouse`: generates plots for experiments performed on mouse synaptosome samples.

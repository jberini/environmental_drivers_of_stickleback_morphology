# Morphology & Environmental Covariates: Ecological Drivers of Stickleback Variation

This repository houses the reproducible analytical pipeline and statistical architectures evaluating how geographic landscape, water chemistry, and zooplankton community gradients drive morphological variation in stickleback populations.

## Repository Directory & Structure

### 💻 Master Scripts
* `Ecological_Drivers_of_Stickleback_Morphological_Variation_GitHUB.Rmd`: The master, fully annotated R Markdown code execution file.
* `Ecological_Drivers_of_Stickleback_Morphological_Variation_GitHUB.md`: The compiled/knitted Markdown document. **Click this file to immediately view all statistical outputs, rendered tables, and figures directly in your browser.**

### 📁 Data Directory (`data/`)
This project utilizes five core relational datasets evaluating phenotypes and ecosystem profiles:
* `stickle_morph.csv`: Master phenotypic data capturing structural and geometric morphometrics across target populations.
* `rf_dat.csv`: Processed environmental metrics configured for Random Forest modeling environments.
* `VancouverSondeData2023.csv`: High-resolution physical water column metrics (DO, pH, temperature, conductivity) captured via data-sonde profiles.
* `ZoopIDDataUpdated.csv`: Taxonomic counts and community composition matrices of local zooplankton prey resources.
* `lake_id_mapping.csv`: Spatial and structural metadata linking independent sample points to geographic lake frameworks.

### 📁 Statistical Summaries & Tables (`figures and tables/`)
Exported comma-separated summary tables of key statistical environments:
* `Table 1. stickleback_regression_table.csv`: Comprehensive summary of univariate linear regression frameworks evaluating environmental drivers against isolated phenotypic traits.
* `Table 2a. dbrda_variance_explained_table.csv`: Distance-based Redundancy Analysis (dbRDA) metrics capturing raw and adjusted cumulative variance ($R^2$) explained by prey communities.
* `Table 2b. dbrda_global_permutation_table.csv`: Global non-parametric permutation test outputs evaluating the full multivariate dbRDA landscape model.
* `Table 2c. dbrda_marginal_axis_table.csv`: Marginal permutation testing evaluating the significance of individual zooplankton PCA axes.

### 📁 Figures & Tables (`figures and tables/`)
High-resolution, publication-ready vector graphics:
* `Figure 1. Stickleback Drivers.pdf`: Main analytical plot mapping primary environmental covariates to morphometric variation.
* `Figure 2. Ecosystem Gradients.pdf`: Ecosystem gradient profiles tracking spatial and raw variation across individual collection sites.

---

## Technical Replication & Requirements

To clone and execute this pipeline locally, open R/RStudio and ensure your local environment contains the following dependencies:

```r
# Required Packages for Execution
install.packages(c("tidyverse", "vegan", "randomForest", "knitr", "ggplot2"))

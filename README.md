# Spatial Epidemiology and Environmental Exposure Modeling using mixed modelling

## Description
This repository contains R scripts for analyzing spatial epidemiological data, environmental exposure, and disease prevalence using mixed-effects models. The analysis integrates exposure data, disease prevalence across ZIP codes, and spatial clustering techniques to assess environmental factors' impact on health outcomes.

The core methodology involves:
- Cleaning and merging exposure, disease, and demographic datasets
- Conducting spatial clustering and Moran's I analysis for spatial autocorrelation
- Implementing Negative Binomial Generalized Linear Mixed Models (GLMM) to assess risk factors
- Parallel computing for efficient processing across multiple health conditions

## Features
- **Data Preprocessing**: Cleans and merges environmental exposure data with disease prevalence
- **Spatial Analysis**: Incorporates Moran’s I and hierarchical clustering for spatial dependency evaluation
- **Mixed-Effects Modeling**: Uses `glmer.nb` to assess relationships between exposure and health outcomes
- **Parallel Computing**: Implements `foreach` and `doParallel` for efficient execution
- **Automated Reporting**: Generates results for each health condition and saves them as CSV files

## Dependencies
Ensure you have the following R packages installed:
```r
install.packages(c("doParallel", "foreach", "lme4", "MuMIn", "readr", "dplyr", "stringr", "tidyr", "ggplot2", "reshape2", "raster", "sf", "spdep", "zipcodeR", "RVAideMemoire", "geosphere"))
```

## Usage
1. Clone this repository:
   ```sh
   git clone https://github.com/yourusername/SpatialEpi_Modeling.git
   cd SpatialEpi_Modeling
   ```
2. Prepare the input datasets (exposure, disease prevalence, and denominators) in CSV format.
3. Run the main script in R:
   ```r
   source("main_script.R")
   ```
4. Results will be saved in the `results/` directory as CSV files.

## Author
Prem Prashant Chaudhary


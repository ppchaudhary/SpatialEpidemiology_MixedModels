#!/bin/bash

#SBATCH --job-name=R_ME_analysis       # Job name
#SBATCH --output=ME_analysis.out       # Output file
#SBATCH --error=ME_analysis.err        # Error file
#SBATCH --cpus-per-task=56             # Number of CPU cores
#SBATCH --gres=lscratch:300            # Scratch space (temporary storage)
#SBATCH --mem=200G                     # Memory limit
#SBATCH --time=120:00:00                # Time limit (hh:mm:ss)

# Load the required R module
module load R/4.2

# Set a library path
LIB_PATH="$HOME/Rlibs"

# Create the directory if it doesn't exist
mkdir -p $LIB_PATH

# Change to the correct directory
cd /data/chaudharyp2/R_ETU/Mixed_effect_to_run/ME_under_18_chunks_set4/df5

# Print the current working directory for debugging
echo "Current working directory:"
pwd

# Install necessary R packages (in case they are not installed)
Rscript -e "packages <- c('data.table', 'doParallel', 'foreach', 'lme4', 'MuMIn', 'readr', 'dplyr', 'stringr', 'tidyr', 'ggplot2', 'reshape2', 'raster', 'sf', 'spdep', 'zipcodeR', 'RVAideMemoire', 'geosphere', 'future');
             new_packages <- packages[!(packages %in% installed.packages()[,'Package'])];
             if(length(new_packages)) install.packages(new_packages, repos='http://cran.r-project.org', lib='$LIB_PATH')"

# Run the R script and redirect output to a log file
Rscript -e ".libPaths(c('$LIB_PATH', .libPaths())); source('ME_AHQR_under_age_v2.R')" > ME_output.log 2>&1

# Check if the R script ran successfully
if [ $? -ne 0 ]; then
    echo "R script encountered an error."
else
    echo "R script completed successfully."
fi

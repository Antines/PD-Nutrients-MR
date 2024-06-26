# MR Analysis Workflow

This repository contains scripts for performing Mendelian Randomization (MR) analysis . The workflow (MR.R) is divided into three primary parts: preparation of the dataset, data joining and formatting, and clumping variants. Each of these parts is automated through a series of functions to facilitate efficient and accurate MR analysis.

## Workflow Overview

### Part 1: Preparation of the Dataset for MR Analysis

This section automates the preparation of phenotype-specific datasets from unique studies based on the traits of interest. The main steps involved are:

1. **Trait Search and Filtering**:
   - The function searches for specified traits within the phenotype/reported trait column (e.g., “B1,”, “protein intake”, Magnesium).
   - Filters the data to include only European ancestry.
   - Counts the overall sample size (including case, control, and replicated individual sample size if present).
   - Selects the highest sample size for each trait.

2. **Data Joining and Formatting**:
   - Joins data by accession total sample size.
   - Formats the data for MR analysis.
   - Calculates the standard error if it is not provided.

3. **Clumping Variants**:
   - Reduces redundancy in the dataset by "clumping" or grouping closely linked SNPs based on physical proximity and linkage disequilibrium (LD) using the following parameters:
     - Physical distance threshold: 10,000 kilobases (kb)
     - LD threshold (measured by the correlation coefficient r²): 0.001
     - P-value threshold for index SNPs: 1
     - P-value threshold for secondary SNPs.

These automated functions streamline the process of preparing and formatting datasets for MR analysis, ensuring consistency and reliability in the data used for subsequent analysis.

## Repository Contents

- `MR.R`: Contains the 3 function for the whole MR worklof icluding preparation of GWAS datasets, MR, Reverse MR.
- `LAVA.R`: Automated LAVA Analisys followed by the LDSC (for more information - https://github.com/josefin-werme/LAVA) 
- `LDSC.sh`: Automated Genetic correlation LDSC test (based on original software mannual - https://github.com/bulik/ldsc)

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Antines/PD-Nutrients-MR.git
   cd PD-Nutrients-MR



# Title: Genome Information Field Theory (GIFT): Using Physics Field Theory to Enhance GWAS Methodology

This study aims to test and validate the GIFT method through in silico modellingm and comparing the performance of GIFT against traditional GWAS methods. We used a variety of approaches to conduct our analysis:

- Genomic Analysis: We used GIFT to identify significant single nucleotide polymorphisms (SNPs)associated with leaf manganese concentration.
- Statistical Analysis: We calculated significance threshold to determine significant SNPs.
- Visualisation: R packages such as 'tidyverse' and 'htmlwidgets' were used to create visual representations of the data, including Manhattan plots to illustrate SNP associations across the Aradopsis Thaliana genome.

Through these techniques, we aim to assess the efectiveness of the GIFT method and explore its potential advantages over traditional GWAS in identifying genetic associations for complex traits.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [Data Description](#Data description)



## Installation 


###Before running any of the .R, .py, or .sh files, you need to set up Conda environments, install the necessary dependencies, and log into ADA. You can do this by following the instructions in Installing_Conda_and_environments. You need to ensure to have the files in 'input files' saved to the correct location.



## Usage

To run the project, follow these steps:

### Step 1: Prepare and Index VCF File

- **Script:** `prep_and_index.sh`
- **Description:** This script prepares the VCF file by indexing the reference genome, compressing and indexing the VCF file, and filling missing allele frequency (AF) tags.
- **Instructions:** Ensure that you have activated the 'combined_env' conda environment, and then run 'prep_and_index.sh'

### Step 2: Perform Analysis

- **Script:** `analysis.sh`
- **Description:** This script performs additional analysis, including running 'GATK VarientsToTable', querying allele frequencies with 'bcftools', and generating statistics
- **Instructions:** Ensure you have activated the 'combined_env' conda environment, run 'analysis.sh', and the script will generate its output.

### Step 3: Clean FASTA File

- **Script:** `cleaning_fasta.sh`
- **Description:** This script cleans the FASTA file by replacing ambiguous nucleotide codes with valid ones
- **Instructions:** Run the `Toumas_PCA.R` script to generate the PCA.

### Step 4: Extracting Allele Frequencies and Solving Duplicates

- **Scripts:** 'create_AF_with_positions.py' 'solve_dupes.py'
- **Description:** These scripts extract allele frequencies and their corresponding positions, and solve duplicate values, writing to a new output file.
- **Instructions:** Run the 'create_AF_with_positions.py' script, followed by the 'solve_dupes.py' script.

### Step 5: Generate SLiM script

- **Script:** 'create_slim_script.py
- **Description:** This python script generates a SLiM script for modelling SNPs in the Aradopsis Thaliana genome. It processes genome and allele frequency data, introduces a beneficial mutation, and assigns phenotypes based on this mutation
- **Instructions:** Run the 'create_slim_script.py' script to create the SLiM script required for simulating genetic scenarios and assigning phenotypes

### Step 5: Run SLiM script

- **Script:** 'run_slim_script.sh'
- **Description:** This script is designed to run the SLiM simulation. 
- **Instructions:** run the 'run_slim_script.sh' script

### Step 6: Extracting Phenotypes from SLiM Output

- **Script:** `extract_phenotypes.py'
- **Description:** Extracting the phenotypic information from our SLiM output
- **Instructions:** Run the `extract_phenotypes.py` script.


### Step 7: Convert SLiM Output to VCF Format

- **Script:** `convert_to_vcf.py`
- **Description:** Converts the SLiM output into the Variant Call Format (VCF)
- **Instructions:** Run the `convert_to_vcf.py' script


### Step 8: GIFT code and sbatch to run 

- **Scripts:** 'physics_GWAS_OOP.py' and 'sbatch_to_run_GIFT.sh
- **Description:** Runs GIFT on our output
- **Instructions:** run 'physics_GWAS_OOP.py', followed by 'sbatch_to_run_GIFT.sh'
- **Accreditations:** The 'physics_GWAS_OOP.py' script, and the template for the 'sbatch_to_run_GIFT.sh' script were provided by Sian Bray

#### note - the next steps are not done on ADA, and use the input data downloaded to your own machine (noted in input data section)

### Step 9: Generating Manhattan Plots from GWAS results
- **Script:** 'manhattan_plot.R'
- **Description:** Processes GIFT results, calculates statistical thresholds, and generates Manhatten plots for visualising significant SNPs
- **Instructions:** Update the path in the script to where the GIFT input files are stored, then proceed to run 'manhatten_plot.R'
- **Accreditations:** The 'physics_GWAS_OOP.py' script, and the template for the 'sbatch_to_run_GIFT.sh' script were provided by Sian Bray

### Step 10: Generating and analysing SNP Data
- **Script:** 'generate_and_analyse_SNP_data.R'
- **Description:** This R script reads our GIFT results and filters for SNPs above the Bonferroni thresholds, counts the total significant SNPs, and identifies the top 5 SNPs with the highest association with leaf ionome manganese concentrations
- **Instructions:** run the 'generate_and_analyse_SNP_data.R' script
- - **Accreditations:** The 'generate_and_analyse_SNP_data.R' script was provided by Sian Bray

### Step 11:
- **Script:** 'SNP_terminal_commands'
- **Description:** This set of commands is used to identify gene IDs for specific SNPs, prepare the GFF3 file for Aradopsis thaliana, intersect SNPs with gene data, and extact unique gene IDs
- **Instructions:** Run each command one by one in terminal



## Credits

### Contributors
- Grace Savanagh
- Sian Bray


### Supervisor
- Sian Bray



## Data Description

### Input Data

All input data for this project has been provided by Sian Bray. The input data required for this project includes the following:

- **ADA:** 'thaliana_genome.fas': a FASTA file containing reference genome sequences
- **ADA:** 'thaliana.vcf': A VCF which contains Aradopsis Thaliana SNP data
- **ADA:** 'physics_GWAS_OOP.py': a python file containing the code to run GIFT

- **Own machine:** leaf_ionome_Mn55_whole_genome_metrics.csv: a CSV file containing the metrics produced by GIFT for leaf ionome concentrations




### Output Data

The output data generated from this project include:

- Manhattan Plots:Visualisations of GIFT results highlighting significant SNPs
- BED Files for Significant SNPs: A BED file containing significant SNPs filtered based on Bonferroni threshold
- Top SNP Details: A list of the top 5 SNPs with the highest association, based on the -log10(p-value) from GWAS results.
- Gene ID Extraction: a list of unique gene IDs extracted from the intersection of SNPs with gene annotations from the GFF3 file


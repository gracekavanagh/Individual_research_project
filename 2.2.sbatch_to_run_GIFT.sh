#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250g
#SBATCH --time=10:00:00
#SBATCH --job-name=SLiM
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/test/error.txt
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk

# Load Conda
source $HOME/.bash_profile

# Activate Conda environment with R
conda activate ngs_pipe_env_gatk4

# Verify R installation
echo "Checking R installation"
Rscript -e "version"

# Verify R packages
echo "Checking installed R packages"
Rscript -e "installed.packages()"

# Run Python script
echo "Running Python script"
python3 /gpfs01/home/mbygk5/individual_project/test/physics_GWAS_OOP.py \
    -v /gpfs01/home/mbygk5/individual_project/test/vcf_for_downstream_analysis/simulation_output_updated.vcf \
    -f /gpfs01/home/mbygk5/individual_project/test/phenotypes.csv \
    -p Phenotype \
    -o /gpfs01/home/mbygk5/individual_project/test/output.csv

# Run R scripts
echo "Running R scripts"
Rscript /gpfs01/home/mbygk5/individual_project/test/output_absolute_theta.R
Rscript /gpfs01/home/mbygk5/individual_project/test/output_pSNP4.R
Rscript /gpfs01/home/mbygk5/individual_project/test/output_pSNP5.R


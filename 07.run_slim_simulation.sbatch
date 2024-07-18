#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250g
#SBATCH --time=10:00:00
#SBATCH --job-name=SLiM
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/assigning_phenotype/SLiM_simulation.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk

# Load necessary modules and activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slim_env  # Adjust to your Conda environment name

# Run the generated SLiM script
slim /gpfs01/home/mbygk5/individual_project/assigning_phenotype/thaliana.slim

# Deactivate Conda environment
conda deactivate

echo "SLiM simulation completed successfully."

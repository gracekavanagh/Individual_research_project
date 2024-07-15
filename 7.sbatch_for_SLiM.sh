#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250g
#SBATCH --time=10:00:00
#SBATCH --job-name=SLiM
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/assigning_phenotype/ERROR_for_combined_thing.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk


# Activate the Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate slim_env


# Run the generated SLiM script
slim /gpfs01/home/mbygk5/individual_project/assigning_phenotype/thaliana.slim

# Deactivate the Conda environment
conda deactivate

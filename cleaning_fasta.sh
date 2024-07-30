#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH --time=01:00:00
#SBATCH --job-name=clean_fasta
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/cleaning_fasta.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk

echo "Starting job at $(date)"

# Define input and output files
INPUT_FASTA="/gpfs01/home/mbygk5/individual_project/thaliana_SNP_analysis/Thaliana_genome.fas"
OUTPUT_FASTA="/gpfs01/home/mbygk5/individual_project/SLiM/Thaliana_genome_clean.fas"

# Copy the input FASTA file to the output location
cp $INPUT_FASTA $OUTPUT_FASTA

# Replace ambiguous nucleotide codes
echo "Replacing ambiguous nucleotides in FASTA file"
sed -i 's/W/A/g' $OUTPUT_FASTA
sed -i 's/R/A/g' $OUTPUT_FASTA
sed -i 's/S/C/g' $OUTPUT_FASTA
sed -i 's/Y/C/g' $OUTPUT_FASTA
sed -i 's/K/T/g' $OUTPUT_FASTA
sed -i 's/M/A/g' $OUTPUT_FASTA
sed -i 's/N/A/g' $OUTPUT_FASTA

echo "Job completed successfully at $(date)"

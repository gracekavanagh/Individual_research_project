#!/usr/bin/env python3

"""
This script generates a SLiM script for modeling Arabidopsis thaliana genome SNPs
based on provided genome and allele frequency data. Additionally, it introduces
a beneficial mutation in the population to simulate linkage and assigns phenotypes
based on the presence of this beneficial mutation.
"""

import os
import time
import random
import numpy as np
import pandas as pd

# Paths to your data files
genome_path = "/gpfs01/home/mbygk5/individual_project/SLiM/Thaliana_genome_clean.fas"
allele_freqs_path = "/gpfs01/home/mbygk5/individual_project/thaliana_SNP_analysis/thaliana_snp_analysis_allele_frequencies.txt"
slim_script_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/thaliana.slim"
full_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/simulation_output.txt"
population_stats_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/population_stats.txt"
ionome_data_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/master_list.csv"

# Parameters
frequency_threshold = 0.9
mutation_limit = 10
beneficial_mutation_pos = 5000  # Position of the beneficial mutation
selection_coefficient = 0.1  # Selective advantage of the beneficial mutation
phenotype_effect = 0.2  # Additional value for phenotype due to beneficial mutation

# Define the selective trait
trait = 'leaf_ionome_Li7'  # Choose the trait for the phenotype

def parse_genome(fasta_path):
    """Parse the genome from a FASTA file."""
    try:
        genome = {}
        current_chrom = None
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_chrom = line[1:]
                    genome[current_chrom] = []
                else:
                    genome[current_chrom].append(line)
        for chrom in genome:
            genome[chrom] = ''.join(genome[chrom])
        return genome
    except Exception as e:
        print(f"Error parsing genome: {e}")
        return None

def read_allele_freqs(freqs_path):
    """Read allele frequencies from a file."""
    try:
        allele_freqs = []
        with open(freqs_path, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) < 3:
                    continue
                try:
                    chrom, pos, freq = fields[0], int(fields[1]), float(fields[2])
                    if freq >= frequency_threshold:
                        allele_freqs.append((chrom, pos, freq))
                except ValueError:
                    continue
        return allele_freqs[:mutation_limit]
    except Exception as e:
        print(f"Error reading allele frequencies: {e}")
        return None

def load_ionome_data(ionome_path):
    """Load the ionome data from a CSV file."""
    try:
        ionome_data = pd.read_csv(ionome_path)
        return ionome_data
    except Exception as e:
        print(f"Error loading ionome data: {e}")
        return None

def assign_background_phenotypes(ionome_data, population_size, trait):
    """Assign background phenotypes based on ionome data."""
    try:
        phenotype_values = ionome_data[trait].values

        # Check and handle NaN values
        if np.any(np.isnan(phenotype_values)):
            print(f"Warning: NaN values found in trait {trait}. Replacing with mean value.")
            phenotype_values = np.nan_to_num(phenotype_values, nan=np.nanmean(phenotype_values))

        mean_phenotype = np.mean(phenotype_values)
        std_phenotype = np.std(phenotype_values)
        background_phenotypes = np.random.normal(mean_phenotype, std_phenotype, population_size)

        # Handle any NaNs that might still be in the array (e.g., if mean/std are NaN)
        if np.any(np.isnan(background_phenotypes)):
            background_phenotypes = np.nan_to_num(background_phenotypes, nan=0.0)

        return background_phenotypes
    except Exception as e:
        print(f"Error assigning background phenotypes: {e}")
        return None

def create_slim_script(genome, allele_freqs, background_phenotypes, output_path):
    """Create a SLiM script based on genome and allele frequencies."""
    start_time = time.time()
    try:
        with open(output_path, 'w', encoding='utf-8') as f:

            f.write("// SLiM script for modeling Arabidopsis thaliana genome SNPs\n")
            f.write("\n")
            f.write("initialize() {\n")
            f.write("    initializeSLiMModelType(\"nonWF\");\n")
            f.write("    initializeMutationRate(7e-9);\n")
            f.write("    initializeMutationType(\"m1\", 0.5, \"f\", 0.0);\n")
            f.write("    initializeGenomicElementType(\"g1\", m1, 1.0);\n")
            f.write("    initializeMutationType(\"m2\", 0.5, \"f\", 0.1);  // Beneficial mutation type\n")

            # Assuming a single chromosome for simplicity
            chrom_length = sum(len(genome[chrom]) for chrom in genome)
            f.write(f"    initializeGenomicElement(g1, 0, {chrom_length - 1});  // First chromosome\n")
            f.write(f"    initializeGenomicElement(g1, {chrom_length}, {2*chrom_length - 1});  // Second chromosome for diploids\n")

            f.write("    initializeRecombinationRate(4e-8);\n")
            f.write(f"    defineConstant(\"phenotype_effect\", {phenotype_effect});\n")
            f.write("}\n\n")

            f.write("1 early() {\n")
            f.write("    sim.addSubpop(\"p1\", 1135);\n")

            for i, phenotype in enumerate(background_phenotypes):
                f.write(f"    p1.individuals[{i}].setValue(\"phenotype\", {phenotype});\n")

            f.write("}\n\n")

            f.write("2 late() {\n")
            f.write("    inds = sim.subpopulations[0].individuals;\n")

            for chrom, pos, freq in allele_freqs:
                pos -= 1
                num_mutants = int(freq * 1135)
                if num_mutants > 0:
                    permuted_indices = random.sample(range(1135), num_mutants)
                    for ind_index in permuted_indices:
                        f.write(f"    inds[{ind_index}].genomes[0].addNewDrawnMutation(m1, {pos}, NULL, NULL);\n")
                        f.write(f"    inds[{ind_index}].genomes[1].addNewDrawnMutation(m1, {pos}, NULL, NULL);\n")

            # Introduce the beneficial mutation and save its object
            f.write("    // Introduce the beneficial mutation\n")
            f.write(f"    beneficial_mutation = inds[0].genomes[0].addNewDrawnMutation(m2, {beneficial_mutation_pos - 1}, NULL, NULL);\n")
            f.write(f"    inds[0].genomes[1].addNewDrawnMutation(m2, {beneficial_mutation_pos - 1}, NULL, NULL);\n")

            f.write("    // Update phenotypes based on beneficial mutation presence\n")
            f.write("    for (ind in inds) {\n")
            f.write("        if (ind.genomes[0].containsMutations(beneficial_mutation)) {\n")
            f.write(f"            ind.setValue(\"phenotype\", ind.getValue(\"phenotype\") + phenotype_effect);\n")
            f.write("        } else if (ind.genomes[1].containsMutations(beneficial_mutation)) {\n")
            f.write(f"            ind.setValue(\"phenotype\", ind.getValue(\"phenotype\") + phenotype_effect);\n")
            f.write("        }\n")
            f.write("    }\n")

            f.write("}\n\n")

            f.write("3 late() {\n")
            f.write(f"    sim.outputFull(\"{full_output_path}\");\n")
            f.write(f"    writeFile(\"{population_stats_path}\", asString(sim.subpopulations[0].individualCount));\n")
            f.write("    sim.simulationFinished();\n")
            f.write("}\n")

        print(f"SLiM script has been created at {output_path}")
        print(f"Time taken: {time.time() - start_time:.2f} seconds")
    except Exception as e:
        print(f"Error creating SLiM script: {e}")

if __name__ == "__main__":
    try:
        genome = parse_genome(genome_path)
        allele_freqs = read_allele_freqs(allele_freqs_path)
        ionome_data = load_ionome_data(ionome_data_path)
        if genome and allele_freqs and ionome_data is not None:
            population_size = 1135
            background_phenotypes = assign_background_phenotypes(ionome_data, population_size, trait)
            if background_phenotypes is not None:
                create_slim_script(genome, allele_freqs, background_phenotypes, slim_script_path)
        else:
            print("Error: Failed to load necessary data.")
    except KeyboardInterrupt:
        print("\nProcess interrupted.")
    except Exception as e:
        print(f"An error occurred: {e}")

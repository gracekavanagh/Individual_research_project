#!/usr/bin/env python3

#generates a SLiM script to model the A.thaliana genome SNPs based on the genome and allele frequency.
#we also introduce a beneficial muttaion and assign a phenotype

import os
import time
import random
import numpy as np
import pandas as pd

#data paths
genome_path = "/gpfs01/home/mbygk5/individual_project/test/Thaliana_genome_clean.fas"
allele_freqs_path = "/gpfs01/home/mbygk5/individual_project/test/allele_freqs_cleaned.txt"
slim_script_path = "/gpfs01/home/mbygk5/individual_project/test/thaliana.slim"
full_output_path = "/gpfs01/home/mbygk5/individual_project/test/simulation_output.txt"
population_stats_path = "/gpfs01/home/mbygk5/individual_project/test/population_stats.txt"
ionome_data_path = "/gpfs01/home/mbygk5/individual_project/test/master_list.csv"


#chosen parameters
frequency_threshold = 0.9
mutation_limit = 10
beneficial_mutation_pos = 5000  #the position of the beneficial mutation
selection_coefficient = 0.1  #the selective advantage of the beneficial mutation
phenotype_effect = 0.2  #the addiitonal value for phenotype due to beneficial mutation

#choosing the benefiical trait
trait = 'leaf_ionome_Li7'  

#parse the genome from FASTA file
def parse_genome(fasta_path):
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

#reading allele frequencies
def read_allele_freqs(freqs_path):
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

#loading the ionome data
def load_ionome_data(ionome_path):
    try:
        ionome_data = pd.read_csv(ionome_path)
        return ionome_data
    except Exception as e:
        print(f"Error loading ionome data: {e}")
        return None

#assigning the background phenotypes based on our ionome data
def assign_background_phenotypes(ionome_data, population_size, trait):
    try:
        phenotype_values = ionome_data[trait].values

        #NAN value handing
        if np.any(np.isnan(phenotype_values)):
            print(f"Warning: NaN values found in trait {trait}. Replacing with mean value.")
            phenotype_values = np.nan_to_num(phenotype_values, nan=np.nanmean(phenotype_values))

        mean_phenotype = np.mean(phenotype_values)
        std_phenotype = np.std(phenotype_values)
        background_phenotypes = np.random.normal(mean_phenotype, std_phenotype, population_size)

        #further handling if mean/sd is still NAN
        if np.any(np.isnan(background_phenotypes)):
            background_phenotypes = np.nan_to_num(background_phenotypes, nan=0.0)

        return background_phenotypes
    except Exception as e:
        print(f"Error assigning background phenotypes: {e}")
        return None

#creating our SLiM script
def create_slim_script(genome, allele_freqs, background_phenotypes, output_path):
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

            #assuming a single chromosome
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

            #introducing the beneficial mutation
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

def main():
    try:
        #loading our genome data
        print("Loading genome data...")
        genome = parse_genome(genome_path)

        #reading allele frequency
        print("Reading allele frequencies...")
        allele_freqs = read_allele_freqs(allele_freqs_path)

        #loading ionome data
        print("Loading ionome data...")
        ionome_data = load_ionome_data(ionome_data_path)

        #assigning background phenotypes
        print("Assigning background phenotypes...")
        background_phenotypes = assign_background_phenotypes(ionome_data, 1135, trait)

        #creating slim script
        print("Creating SLiM script...")
        create_slim_script(genome, allele_freqs, background_phenotypes, slim_script_path)

        print("SLiM script generation completed successfully.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()

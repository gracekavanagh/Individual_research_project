#!/usr/bin/env python3

import pandas as pd

# paths to SLiM output and new phenotype file
simulation_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/simulation_output.txt"
phenotypes_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/phenotypes.csv"

#extracting the phenotypic data
def extract_phenotypes(simulation_output_path, phenotypes_output_path):
    phenotypes = []

    with open(simulation_output_path, 'r') as f:
        for line in f:
            # Assuming the output format of SLiM stores phenotype values in some identifiable format
            # For this example, let's assume phenotypes are marked with a specific prefix like "PHEN:"
            if line.startswith("PHEN:"):
                fields = line.strip().split()
                try:
                    ind_id = int(fields[1])  # assuming second field is Individual ID
                    phenotype = float(fields[2])  # assuming third field is phenotype value
                    phenotypes.append((ind_id, phenotype))
                except ValueError:
                    continue

    phenotypes_df = pd.DataFrame(phenotypes, columns=['Individual_ID', 'Phenotype'])
    phenotypes_df.to_csv(phenotypes_output_path, index=False)

def main():
    extract_phenotypes(simulation_output_path, phenotypes_output_path)

if __name__ == "__main__":
    main()





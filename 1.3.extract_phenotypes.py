#!/usr/bin/env python3

"""
This script extracts phenotype information from the SLiM simulation output
and saves it into a CSV file.
"""

import pandas as pd

# Paths to the SLiM simulation output and the phenotype CSV output
simulation_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/simulation_out$
phenotypes_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/phenotypes.csv"

def extract_phenotypes(simulation_output_path, phenotypes_output_path):
    """Extract phenotypes from the SLiM simulation output and save to a CSV file."""
    phenotypes = []

    with open(simulation_output_path, 'r') as f:
        for line in f:
            if line.startswith("#OUT:"):
                continue
            fields = line.strip().split()
            if len(fields) >= 2:
                try:
                    ind_id = int(fields[0])
                    phenotype = float(fields[1])
                    phenotypes.append((ind_id, phenotype))
                except ValueError:
                    continue

    phenotypes_df = pd.DataFrame(phenotypes, columns=['Individual_ID', 'Phenotype'])
    phenotypes_df.to_csv(phenotypes_output_path, index=False)
    print(f"Phenotypes saved to {phenotypes_output_path}")

def main():
    extract_phenotypes(simulation_output_path, phenotypes_output_path)

if __name__ == "__main__":
    main()




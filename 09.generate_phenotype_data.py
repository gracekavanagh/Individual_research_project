#need to have done: pip install pandas

import pandas as pd
import re

# Path to your SLiM output file
slim_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/simulation_output.txt"
phenotype_output_path = "/gpfs01/home/mbygk5/individual_project/step_1_GIFT_and_GWAS/phenotypes.txt"

# Read SLiM output file
with open(slim_output_path, 'r') as f:
    lines = f.readlines()

# Function to check if a string represents a valid integer
def is_valid_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

# Parse phenotype data
phenotype_data = []
for line in lines:
    if line.startswith("p1:"):  # Check if line starts with "p1:"
        parts = line.strip().split()
        id_part = parts[0].split(":")[1]
        if is_valid_int(id_part) and len(parts) > 2:  # Check if the ID part is a valid integer and there are enough parts
            try:
                ind_id = int(id_part)  # Extract numeric ID after "p1:"
                phenotypes = list(map(float, parts[2:]))  # Convert remaining parts to float as phenotypes
                phenotype_data.append((ind_id, *phenotypes))  # Append ID and phenotypes as separate columns
            except (ValueError, IndexError) as e:
                print(f"Skipping line due to error: {line.strip()}\nError: {e}")
        else:
            print(f"Skipping line due to unexpected format: {line.strip()}")

# Check the parsed data structure
print("Parsed phenotype data (first 5 entries):", phenotype_data[:5])  # Print first 5 entries for verification

# Convert to DataFrame if there is valid data
if phenotype_data:
    phenotype_df = pd.DataFrame(phenotype_data)
    # Add column names dynamically based on number of phenotype columns
    columns = ["IID"] + [f"Phenotype_{i+1}" for i in range(phenotype_df.shape[1] - 1)]
    phenotype_df.columns = columns

    # Add family ID (FID) column
    phenotype_df["FID"] = 1  # Assuming all individuals belong to the same family

    # Reorder columns to FID, IID, and Phenotypes
    cols_order = ["FID", "IID"] + [col for col in phenotype_df.columns if col.startswith("Phenotype")]
    phenotype_df = phenotype_df[cols_order]

    # Save to a phenotype file
    phenotype_df.to_csv(phenotype_output_path, sep="\t", index=False)
    print(f"Phenotype data saved to {phenotype_output_path}")
else:
    print("No valid phenotype data found.")

import pandas as pd

# Load the phenotype data
phenotype_df = pd.read_csv("/gpfs01/home/mbygk5/individual_project/assigning_phenotype/phenotypes.csv")

# Display basic information
print(phenotype_df.head())  # Display the first few rows
print(phenotype_df.describe())  # Display summary statistics
print(phenotype_df.info())  # Display DataFrame info

# Check for any missing values
missing_values = phenotype_df.isnull().sum()
print(f"Missing values:\n{missing_values}")

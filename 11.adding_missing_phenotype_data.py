import pandas as pd
import numpy as np

# Paths to your data files
phenotype_data_path = '/gpfs01/home/mbygk5/individual_project/step_1_GIFT_and_GWAS/phenotypes.txt'
reference_allele_data_path = '/gpfs01/home/mbygk5/individual_project/test/table.txt'

# Load phenotype data
phenotype_data = pd.read_csv(phenotype_data_path, sep='\t')

# Check column names
print("Columns in phenotype data:", phenotype_data.columns)

# Load reference allele data
reference_allele_data = pd.read_csv(reference_allele_data_path, sep=r'\s+')

# Identify reference alleles
reference_allele_data['Reference_Allele'] = reference_allele_data['AF'] < 0.5
reference_allele_mapping = reference_allele_data.set_index('POS')['Reference_Allele'].to_dict()

# Merge data
if 'IID' not in phenotype_data.columns:
    print("Error: 'IID' column not found in phenotype data.")
    raise KeyError("'IID' column is missing from the phenotype data.")

iid_list = phenotype_data['IID'].unique()
reference_allele_df = pd.DataFrame({'IID': iid_list})
reference_allele_df['Reference_Allele'] = reference_allele_df['IID'].apply(lambda x: reference_allele_mapping.get(x, False))
merged_data = pd.merge(phenotype_data, reference_allele_df, on='IID', how='left')

# Calculate and impute missing values
def calculate_reference_mean(phenotype_column, merged_data):
    reference_individuals = merged_data['Reference_Allele'] == True
    reference_phenotypes = merged_data.loc[reference_individuals, phenotype_column]
    return reference_phenotypes.mean()

for column in phenotype_data.columns:
    if column not in ['FID', 'IID'] and phenotype_data[column].isna().sum() > 0:
        reference_mean_value = calculate_reference_mean(column, merged_data)
        if pd.isna(reference_mean_value):
            print(f"Warning: Reference mean value for {column} is NaN. Skipping imputation for this column.")
        else:
            print(f"Imputing missing values in {column} with reference allele mean value {reference_mean_value}")
            phenotype_data[column] = phenotype_data[column].fillna(reference_mean_value)

# Save cleaned data
phenotype_txt_path = '/gpfs01/home/mbygk5/individual_project/step_1_GIFT_and_GWAS/complete_phenotype.txt'
phenotype_data.to_csv(phenotype_txt_path, sep='\t', index=False)
print(f"Cleaned phenotype data saved to {phenotype_txt_path}")

# Define input and output file paths
input_vcf = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/vcf_for_downstream_analysis$
output_vcf = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/vcf_for_downstream_analysi$

# Open input and output files
with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
    for line in fin:
        line = line.strip()
        if line.startswith('#CHROM'):
            # Modify the header line to add FORMAT and sample column
            header_parts = line.split('\t')
            header_parts.extend(['FORMAT', 'Sample1'])  # Replace Sample1 with your actual sample n$
            fout.write('\t'.join(header_parts) + '\n')
        elif not line.startswith('#'):
            # Process variant lines
            parts = line.split('\t')
            # Append '.' for FORMAT and GT field under Sample1 column
            parts.extend(['.', '.'])
            fout.write('\t'.join(parts) + '\n')
        else:
            # Write other header lines as they are
            fout.write(line + '\n')

print(f"Updated VCF written to {output_vcf}")





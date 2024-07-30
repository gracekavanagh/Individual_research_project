input_table = "/gpfs01/home/mbygk5/individual_project/test/table.txt"
output_freqs = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/allele_freqs_with_positions.txt"

with open(input_table, 'r') as infile, open(output_freqs, 'w') as outfile:
    header = infile.readline()  # Skip the header
    for line in infile:
        parts = line.strip().split()
        if len(parts) < 8:
            continue  # Skip lines that don't have enough columns
        pos = int(parts[1])
        af = float(parts[7])  # Adjust index if necessary
        outfile.write(f"{pos} {af}\n")

print(f"Processed allele frequencies to {output_freqs}")

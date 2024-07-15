# Path to your simulation output file
simulation_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/simulation_output.txt"
# Path to your VCF output file
vcf_output_path = "/gpfs01/home/mbygk5/individual_project/assigning_phenotype/vcf_for_downstream_analysis/simulation_output.vcf"

def parse_simulation_output(simulation_output_path):
    """Parse simulation output file to extract mutation information."""
    mutations = []
    in_mutations_section = False
    with open(simulation_output_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Mutations:"):
                in_mutations_section = True
                continue
            if in_mutations_section:
                if line == "" or line.startswith("Individuals:"):
                    break  # End of mutations section or start of individual data
                parts = line.split()
                try:
                    mut_id = parts[0]
                    chrom = "1"  # Default chromosome, adjust if needed
                    pos = parts[3]  # Assuming the position is the 4th element
                    ref = "A"  # Placeholder reference allele
                    alt = "T"  # Placeholder alternate allele
                    mutations.append((chrom, pos, mut_id, ref, alt))
                except IndexError:
                    print(f"Error parsing line: {line}")
                    continue
    return mutations

def write_vcf(mutations, vcf_output_path):
    """Write mutations to a VCF file."""
    with open(vcf_output_path, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write mutation data
        for chrom, pos, mut_id, ref, alt in mutations:
            f.write(f"{chrom}\t{pos}\t{mut_id}\t{ref}\t{alt}\t.\t.\t.\n")

# Parse simulation output and write to VCF
mutations = parse_simulation_output(simulation_output_path)
write_vcf(mutations, vcf_output_path)

print(f"VCF file has been created at {vcf_output_path}")

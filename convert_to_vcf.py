import csv
import random

#selecting our paths
slim_output_path = "/gpfs01/home/mbygk5/individual_project/test/simulation_output.txt"
vcf_output_path = "/gpfs01/home/mbygk5/individual_project/test/simulation_output.vcf"

# reading our SLiM output
def read_slim_output(slim_output_path):
    individuals = []
    with open(slim_output_path, 'r') as file:
        for line in file:
            if line.startswith("IND"):
                data = line.strip().split()
                individuals.append({
                    "id": data[1],
                    "genotype": data[2:],  
                    "phenotype": data[-1]  
                })
    return individuals

# creating VCF
def create_vcf(individuals):
    vcf_content = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([ind['id'] for ind in individuals])
    ]
    for i, ind in enumerate(individuals):
        vcf_content.append(f"1\t{100+i}\t.\tA\tT\t.\t.\t.\tGT\t" + "\t".join(["0/1" for _ in individuals]))
    return "\n".join(vcf_content)

def write_vcf(vcf_output_path, vcf_content):
    with open(vcf_output_path, 'w') as file:
        file.write(vcf_content)


individuals = read_slim_output(slim_output_path)
vcf_content = create_vcf(individuals)
write_vcf(vcf_output_path, vcf_content)




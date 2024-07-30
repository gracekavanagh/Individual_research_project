
#############################################run this to prep vcf 

#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=10:00:00
#SBATCH --job-name=thaliana_SNP_analysis
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/test/thaliana_SNP_analysis.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk

echo "Starting job at $(date)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate combined_env

# Paths
REFERENCE="/gpfs01/home/mbygk5/individual_project/test/Thaliana_genome.fas"
VCF="/gpfs01/home/mbygk5/individual_project/test/thaliana.vcf"
COMPRESSED_VCF="${VCF}.gz"
INDEXED_VCF="${COMPRESSED_VCF}.tbi"
OUTPUT_TABLE="/gpfs01/home/mbygk5/individual_project/test/table.txt"


echo "Starting job at $(date)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate combined_env

# Paths
REFERENCE="/gpfs01/home/mbygk5/individual_project/test/Thaliana_genome.fas"
VCF="/gpfs01/home/mbygk5/individual_project/test/thaliana.vcf"
COMPRESSED_VCF="${VCF}.gz"
INDEXED_VCF="${COMPRESSED_VCF}.tbi"
OUTPUT_TABLE="/gpfs01/home/mbygk5/individual_project/test/table.txt"

# Step 1: Ensure the reference genome is indexed
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Indexing reference genome"
    samtools faidx $REFERENCE
fi

# Step 2: Compress the VCF file using bgzip
if [ ! -f "$COMPRESSED_VCF" ]; then
    echo "Compressing VCF file with BGZF"
    bgzip -c $VCF > $COMPRESSED_VCF
fi

# Step 3: Index the compressed VCF file using tabix
if [ ! -f "$INDEXED_VCF" ]; then
    echo "Indexing VCF file"
    tabix -p vcf $COMPRESSED_VCF
fi

# Step 4: Fill missing AF tags in VCF file using bcftools
if [ ! -f "filled_af.vcf.gz" ]; then
    echo "Filling missing AF tags in VCF file"
    bcftools +fill-tags $COMPRESSED_VCF -Oz -o filled_af.vcf.gz -- -t AF
    tabix -p vcf filled_af.vcf.gz
fi

echo "Job completed successfully at $(date)"

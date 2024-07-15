

#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=10:00:00
#SBATCH --job-name=thaliana_SNP_analysis
#SBATCH --output=/gpfs01/home/mbygk5/individual_project/thaliana_SNP_analysis/thaliana_SNP_analysis.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbygk5@exmail.nottingham.ac.uk

echo "Starting job at $(date)"

# Load necessary modules
source /gpfs01/home/mbygk5/miniconda3/bin/activate
echo "Activating combined_env"
conda activate combined_env

if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment"
    exit 1
fi

# Define input and output files
VCF_FILE='/gpfs01/home/mbygk5/individual_project/thaliana_SNP_analysis/filled_af.vcf.gz'
REF_GENOME='Thaliana_genome.fas'
OUTPUT_PREFIX='/gpfs01/home/mbygk5/individual_project/thaliana_SNP_analysis/thaliana_snp_analysis'

# Index the reference genome if not already indexed
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Indexing reference genome"
    samtools faidx $REF_GENOME
    if [ $? -ne 0 ]; then
        echo "Error indexing reference genome"
        exit 1
    fi
else
    echo "Reference genome already indexed"
fi

# Check for existing VCF index files
if [ ! -f "${VCF_FILE}.csi" ] && [ ! -f "${VCF_FILE}.tbi" ]; then
    echo "Indexing VCF file"
    bcftools index -t $VCF_FILE  # Create a .tbi index
    if [ $? -ne 0 ]; then
        echo "Error indexing VCF file"
        exit 1
    fi
else
    echo "VCF file already indexed"
fi


# Running GATK VariantsToTable with debug options
echo "Running GATK VariantsToTable"
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantsToTable \
    -V $VCF_FILE \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F DP \
    -O ${OUTPUT_PREFIX}_variants_table.txt
if [ $? -ne 0 ]; then
    echo "Error running GATK VariantsToTable"
    exit 1
fi

# Running BCFtools query
echo "Running BCFtools query"
bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF_FILE > ${OUTPUT_PREFIX}_allele_frequencies.txt
if [ $? -ne 0 ]; then
    echo "Error running BCFtools query"
    exit 1
fi

# Running BCFtools stats
echo "Running BCFtools stats"
bcftools stats $VCF_FILE > ${OUTPUT_PREFIX}_bcftools_stats.txt
if [ $? -ne 0 ]; then
    echo "Error running BCFtools stats"
    exit 1
fi

echo "Job completed successfully at $(date)"

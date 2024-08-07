

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



# input and output files
VCF_FILE='/gpfs01/home/mbygk5/individual_project/test/filled_af.vcf.gz'
REF_GENOME='/gpfs01/home/mbygk5/individual_project/test/Thaliana_genome.fas'
OUTPUT_PREFIX='/gpfs01/home/mbygk5/individual_project/test/thaliana_snp_analysis'

# index reference genome
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Indexing reference genome"
    samtools faidx $REF_GENOME
fi

# index VCF
if [ ! -f "${VCF_FILE}.tbi" ]; then
    echo "Indexing VCF file"
    bcftools index -t $VCF_FILE
fi

# VarientsToTables to make our data easier to read
echo "Running GATK VariantsToTable"
gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantsToTable \
    -V $VCF_FILE \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F DP \
    -O ${OUTPUT_PREFIX}_variants_table.txt

#  BCFtools query and stats
bcftools query -f '%CHROM\t%POS\t%AF\n' $VCF_FILE > ${OUTPUT_PREFIX}_allele_frequencies.txt
bcftools stats $VCF_FILE > ${OUTPUT_PREFIX}_bcftools_stats.txt



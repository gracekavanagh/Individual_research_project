#R code#


#generate list of genes from GIFT focusing on pSNP5 ----

# Load necessary libraries
library("tidyverse")
library("htmlwidgets")

#read GWAS results
GWAS_result1 <- read.csv("/Users/gracekavanagh/Downloads/leaf_ionome_data/leaf_ionome_Mn55_whole_genome_metrics.csv")
GWAS_result1[GWAS_result1 == Inf] <- NA
GWAS_result1[GWAS_result1 == -Inf] <- NA
GWAS_result <- GWAS_result1[complete.cases(GWAS_result1), ]

#calculate the Bonferroni threshold
bt <- 0.01 / (nrow(GWAS_result) * 1135) # times max number of tests per p-value
bf_thres <- -log10(bt)

#filter SNPs above the Bonferroni threshold
significant_snps <- GWAS_result %>% filter(-log10(pSNP5) > bf_thres)

#check the column names to ensure we are using the correct ones
colnames(significant_snps)

#convert to BED format
bed_data <- significant_snps %>%
  select(CHROM, POS) %>%
  mutate(start = POS, end = POS + 1) %>%
  select(CHROM, start, end)

#write out to BED file
write.table(bed_data, file = "/Users/gracekavanagh/Downloads/leaf_ionome_data/significant_snps.bed", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# significant SNPs ----

#count the number of significant SNPs ----
num_significant_snps <- nrow(significant_snps)
num_significant_snps

#count the number of significant SNPs per chromosome
significant_snps_per_chromosome <- significant_snps %>%
  group_by(CHROM) %>%
  summarize(count = n())

#print the counts
print(significant_snps_per_chromosome)




#finding top 5 SNPs from manhatten with highest association ----

#load the necessary library
library(dplyr)

#read and clean the GWAS results
GWAS_result <- read.csv("/Users/gracekavanagh/Downloads/leaf_ionome_data/leaf_ionome_Mn55_whole_genome_metrics.csv")
GWAS_result[GWAS_result == Inf] <- NA
GWAS_result[GWAS_result == -Inf] <- NA
GWAS_result <- GWAS_result[complete.cases(GWAS_result), ]

#add a column for -log10(p-value) if not already present
GWAS_result <- GWAS_result %>%
  mutate(log_p_value = -log10(pSNP5))

#identify the top 5 SNPs with the highest -log10(p-value)
top_snps <- GWAS_result %>%
  arrange(desc(log_p_value)) %>%
  slice(1:5)

#print the details of the top 5 SNPs
print(top_snps)



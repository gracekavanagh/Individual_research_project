
#paths to files
allele_freqs_file = "/gpfs01/home/mbygk5/individual_project/test/allele_freqs_with_positions.txt"
cleaned_file = "/gpfs01/home/mbygk5/individual_project/test/allele_freqs_cleaned.txt"

freq_dict = defaultdict(list)

#reading the AF file and aggregating frequencies
with open(allele_freqs_file, 'r') as infile:
    for line in infile:
        parts = line.strip().split()
        if len(parts) == 2:
            pos, freq = int(parts[0]), float(parts[1])
            freq_dict[pos].append(freq)

#cleaning frequencies
with open(cleaned_file, 'w') as outfile:
    for pos, freqs in sorted(freq_dict.items()):
        avg_freq = sum(freqs) / len(freqs)  # Average the frequencies
        outfile.write(f"{pos} {avg_freq:.8f}\n")



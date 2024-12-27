import random

def select_subset_of_snps(input_file, output_file, subset_size=100000):
    with open(input_file, 'r') as infile:
        snps = [line.strip() for line in infile.readlines()]
    
    subset_size = min(subset_size, len(snps))  # Adjust subset size if fewer SNPs are available
    selected_snps = random.sample(snps, subset_size)
    
    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(selected_snps) + '\n')
    
    print(f"Selected {subset_size} SNPs and saved to {output_file}")

def filter_vcf_by_positions(vcf_file, pos_file, output_vcf):
    # Read positions from the position text file
    with open(pos_file, 'r') as pos_file:
        positions_to_include = set(line.strip() for line in pos_file.readlines())

    # Open the VCF file and prepare to write the filtered output
    with open(vcf_file, 'r') as vcf:
        with open(output_vcf, 'w') as output:
            # Write the header lines (lines that start with ## or #)
            for line in vcf:
                if line.startswith('##') or line.startswith('#'):
                    output.write(line)
                else:
                    # For the data rows, check if the position is in the list of positions to include
                    fields = line.split('\t')
                    position = fields[1]  # VCF position is in the second column
                    if position in positions_to_include:
                        output.write(line)
    
    print(f"Filtered VCF saved to {output_vcf}")


import csv
import random

def filter_vcf_by_samples(input_tsv, input_vcf, output_vcf, output_tsv, n_samples):
    """
    Filters a VCF file based on a list of randomly selected sample names provided in a TSV file.
    Only the samples listed in the TSV file are included in the output VCF.
    
    Randomly selects 'n_samples' samples from the TSV file and writes their IDs and classes
    to a new TSV file.

    Args:
        input_tsv (str): Path to the input TSV file containing the sample names and their classes (two columns).
        input_vcf (str): Path to the input VCF file to be filtered.
        output_vcf (str): Path to the output VCF file containing only the selected samples.
        output_tsv (str): Path to the output TSV file containing the selected sample IDs and their classes.
        n_samples (int): Number of samples to randomly select from the TSV file.

    Returns:
        None. The function writes the filtered VCF to the specified output VCF file and the selected samples to the output TSV file.
    
    Example:
        filter_vcf_by_samples('samples_to_include.tsv', 'input_data.vcf', 'filtered_data.vcf', 'filtered_samples.tsv', 50)
    """
    
    # Read the sample names and their corresponding classes from the TSV file
    samples_dict = {}
    with open(input_tsv, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            sample_id, sample_class = row[0], row[1]
            samples_dict[sample_id] = sample_class
    
    # Randomly select 'n_samples' from the list of sample names
    if n_samples > len(samples_dict):
        raise ValueError(f"Requested n_samples ({n_samples}) is greater than the available samples ({len(samples_dict)})")
    
    selected_samples = random.sample(list(samples_dict.keys()), n_samples)
    
    # Create a new TSV file with selected sample IDs and their classes
    with open(output_tsv, 'w', newline='') as out_tsv:
        writer = csv.writer(out_tsv, delimiter='\t')
        for sample_id in selected_samples:
            writer.writerow([sample_id, samples_dict[sample_id]])
    
    # Open the input VCF and process the file
    with open(input_vcf, 'r') as vcf_file, open(output_vcf, 'w') as out_vcf:
        header_written = False
        for line in vcf_file:
            # Handle header lines (starting with '#')
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # Extract the sample names from the header
                    header_fields = line.strip().split('\t')
                    all_samples = header_fields[9:]  # Sample names start from index 9
                    # Filter the header to include only selected samples
                    filtered_samples = header_fields[:9] + [sample for sample in all_samples if sample in selected_samples]
                    out_vcf.write('\t'.join(filtered_samples) + '\n')
                else:
                    # For other header lines, just write them as-is
                    out_vcf.write(line)
                continue
            
            # Process data lines (not header lines)
            fields = line.strip().split('\t')
            sample_names = fields[9:]
            
            # Filter the genotypes to include only those for the selected samples
            filtered_genotypes = fields[:9]  # Keep the first 9 columns (standard INFO columns)
            for i, sample in enumerate(sample_names):
                if all_samples[i] in selected_samples:
                    filtered_genotypes.append(sample_names[i])

            # Write the filtered SNP information to the output VCF
            out_vcf.write('\t'.join(filtered_genotypes) + '\n')

    print(f"Filtered VCF file written to {output_vcf}")
    print(f"Filtered TSV file written to {output_tsv}")



# Example usage
# select_subset_of_snps(r'E:\GeMorph\Ancestry\SALAI-Net\data\ref.txt', 'data/final/selected_snps.txt')

# filter_vcf_by_positions(r'E:\GeMorph\Ancestry\SALAI-Net\data\final\ref_panel.vcf', r'E:\GeMorph\Ancestry\SALAI-Net\data\final\selected_snps.txt', 'data/final/ref_panel_chr22.vcf')

# filter_vcf_by_positions(r'E:\GeMorph\Ancestry\SALAI-Net\data\final\train_dataset.vcf', r'E:\GeMorph\Ancestry\SALAI-Net\data\final\selected_snps.txt', 'data/final/train.vcf')

filter_vcf_by_positions(r'E:\GeMorph\Ancestry\SALAI-Net\data\final\test_dataset.vcf', r'E:\GeMorph\Ancestry\SALAI-Net\data\final\selected_snps.txt', 'data/final/test.vcf')

# filter_vcf_by_samples(r'E:\GeMorph\Ancestry\SALAI-Net\data\final\ref_panel_map.tsv', r'E:\GeMorph\Ancestry\SALAI-Net\data\final\ref_panel_chr22.vcf', 'small_ref_panel_chr22.vcf', 'small_ref_panel.tsv', 100)

import csv
import random

def extract_samples_and_header(input_vcf, output_vcf, output_tsv, input_tsv):
    # Initialize variables
    selected_samples = []
    sample_classes = {}

    # Read the sample classes from the input TSV file
    with open(input_tsv, 'r') as tsv_file:
        print("reading input map")
        reader = csv.reader(tsv_file, delimiter='\t')
        # Skip the header row if there is one (use next(reader))
        for row in reader:
            sample_name, class_label = row
            if class_label not in sample_classes:
                sample_classes[class_label] = []
            sample_classes[class_label].append(sample_name)

    # Ensure each class has exactly 100 samples
    selected_sample_names = []
    for class_label, samples in sample_classes.items():
        if len(samples) > 100:
            # Randomly select 100 samples
            selected_samples_for_class = random.sample(samples, 100)
        else:
            # If there are fewer than 100 samples, select all
            selected_samples_for_class = samples
        selected_sample_names.extend(selected_samples_for_class)
    
    # Open input VCF file and output VCF file
    with open(input_vcf, 'r') as vcf_file, open(output_vcf, 'w') as out_vcf:
        print("reading input vcf")
        
        # Process the file line by line
        header_lines = []
        all_sample_names = []
        
        # Process the header lines and extract sample names
        for line in vcf_file:
            if line.startswith('#'):
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    # The sample names are in the last header line starting with '#CHROM'
                    all_sample_names = line.strip().split('\t')[9:]
            else:
                break  # Data section starts here
        
        # Write all header lines to the output VCF file
        out_vcf.writelines(header_lines)
        
        # Process the VCF data section line by line
        for line in vcf_file:
            fields = line.strip().split('\t')
            selected_columns = fields[:9]  # Keep the first 9 columns (metadata columns)
            
            # Add data for selected samples only
            for i, sample_name in enumerate(all_sample_names):
                if sample_name in selected_sample_names:
                    selected_columns.append(fields[i + 9])  # Sample data starts from the 10th column
            
            # Write the filtered data line to the output VCF
            out_vcf.write('\t'.join(selected_columns) + '\n')

    print("writing output tsv")
    # Write the selected sample names and their classes to the TSV file
    with open(output_tsv, 'w', newline='') as tsv_out:
        writer = csv.writer(tsv_out, delimiter='\t')

        # Write each sample and its class incrementally
        for sample_name in selected_sample_names:
            # Find the class of the sample
            for class_label, samples in sample_classes.items():
                if sample_name in samples:
                    writer.writerow([sample_name, class_label])
                    break

    print("Ref_panel created")

input_vcf = r'E:\GeMorph\Ancestry\SALAI-Net\data\utils\ch_22_filtered.vcf'
output_vcf = r'ref_panel_chr22.vcf'
output_tsv = 'ref_panel_map.tsv'
input_tsv = r'E:\GeMorph\Ancestry\SALAI-Net\data\1kmap.tsv'
extract_samples_and_header(input_vcf, output_vcf, output_tsv, input_tsv)
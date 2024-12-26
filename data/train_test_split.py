import csv
import random
from typing import Set


def compare_and_filter_tsv(input_tsv1, input_tsv2, output_tsv):
    """
    Compares two TSV files based on the first column (sample names) and filters out the samples
    from the first TSV that are also present in the second TSV. The remaining samples from the first
    TSV are written to a new output TSV file.
    
    Args:
        input_tsv1 (str): Path to the first TSV file (contains sample names and associated labels).
        input_tsv2 (str): Path to the second TSV file (contains sample names to be filtered out).
        output_tsv (str): Path to the output TSV file where the filtered samples from the first TSV will be written.
        
    Returns:
        None. The function writes the filtered results to the output_tsv file.
    
    Example:
        compare_and_filter_tsv('file1.tsv', 'file2.tsv', 'filtered_output.tsv')
        
    In the above example, the function will filter out any rows from 'file1.tsv' whose first column value 
    (sample name) is also found in 'file2.tsv' and write the remaining rows to 'filtered_output.tsv'.
    """
    
    # Read the second TSV file into a set of sample names (first column only)
    with open(input_tsv2, 'r') as tsv_file2:
        samples_in_tsv2 = {row[0] for row in csv.reader(tsv_file2, delimiter='\t')}
    
    # Read the first TSV and filter out the samples that are in the second TSV
    filtered_samples = [
        row for row in csv.reader(open(input_tsv1, 'r'), delimiter='\t')
        if row[0] not in samples_in_tsv2
    ]
    
    # Write filtered samples to output TSV
    with open(output_tsv, 'w', newline='') as tsv_out:
        writer = csv.writer(tsv_out, delimiter='\t')
        writer.writerows(filtered_samples)

    print(f"Filtered samples written to {output_tsv}")


def stratified_split(input_tsv, train_tsv, test_tsv, ratio=0.8):
    """
    Splits samples from the input TSV into stratified train and test sets, ensuring that each class 
    is represented in both sets according to the specified ratio.

    Args:
        input_tsv (str): Path to the input TSV file containing sample names and their class labels.
        train_tsv (str): Path to the output TSV file for the training set.
        test_tsv (str): Path to the output TSV file for the test set.
        ratio (float): Proportion of data to be used for training (default is 0.8).

    Returns:
        None. The function writes the stratified train and test sets to the specified output TSV files.
    
    Example:
        stratified_split('samples.tsv', 'train_samples.tsv', 'test_samples.tsv', ratio=0.8)
    """
    
    # Read the samples and their classes from the input TSV
    class_samples = {}
    with open(input_tsv, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            sample_name, class_label = row
            if class_label not in class_samples:
                class_samples[class_label] = []
            class_samples[class_label].append(row)

    # Stratified Split - ensuring all classes are in both train and test sets
    train_samples = []
    test_samples = []
    
    for class_label, samples in class_samples.items():
        random.shuffle(samples)  # Shuffle the samples within the class
        split_index = int(len(samples) * ratio)
        
        # Split the class samples into train and test sets
        train_samples.extend(samples[:split_index])
        test_samples.extend(samples[split_index:])
    
    # Write the train samples to the train TSV file
    with open(train_tsv, 'w', newline='') as train_file:
        writer = csv.writer(train_file, delimiter='\t')
        writer.writerow(['Sample_Name', 'Class'])  # Write header
        writer.writerows(train_samples)

    # Write the test samples to the test TSV file
    with open(test_tsv, 'w', newline='') as test_file:
        writer = csv.writer(test_file, delimiter='\t')
        writer.writerow(['Sample_Name', 'Class'])  # Write header
        writer.writerows(test_samples)

    print(f"Stratified Train and Test sets written to {train_tsv} and {test_tsv}")


def filter_vcf_by_samples(input_tsv: str, input_vcf: str, output_vcf: str) -> None:
    """
    Filters a VCF file based on a list of sample names provided in a TSV file.
    Uses buffered writing for improved I/O performance.

    Args:
        input_tsv (str): Path to the input TSV file containing sample names.
        input_vcf (str): Path to the input VCF file to be filtered.
        output_vcf (str): Path to the output VCF file.
    """
    # Read sample names into a set for O(1) lookups
    with open(input_tsv, 'r') as tsv_file:
        samples_to_include: Set[str] = {row[0] for row in csv.reader(tsv_file, delimiter='\t')}
    
    # Use a larger buffer size for writing (1MB)
    buffer_size = 1024 * 1024  # 1MB buffer
    
    with open(input_vcf, 'r') as vcf_file, \
         open(output_vcf, 'w', buffering=buffer_size) as out_vcf:
        
        sample_indices = []
        header_written = False
        
        for line in vcf_file:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # Process header line with sample names
                    header_fields = line.strip().split('\t')
                    all_samples = header_fields[9:]
                    
                    # Calculate indices of samples to keep (done once)
                    sample_indices = [i for i, sample in enumerate(all_samples) 
                                   if sample in samples_to_include]
                    
                    # Write filtered header
                    filtered_header = (header_fields[:9] + 
                                    [all_samples[i] for i in sample_indices])
                    out_vcf.write('\t'.join(filtered_header) + '\n')
                else:
                    # Write other header lines as-is
                    out_vcf.write(line)
                continue
            
            # Process data lines using pre-calculated indices
            fields = line.rstrip().split('\t')
            filtered_line = (fields[:9] +  # Fixed columns
                           [fields[i + 9] for i in sample_indices])  # Selected samples
            out_vcf.write('\t'.join(filtered_line) + '\n')

    print(f"Filtered VCF file written to {output_vcf}")


#Usage

# get the samples which are not in the ref
compare_and_filter_tsv('1kmap.tsv', 'test.tsv', 'test_tsv.tsv')
# divide the samples into train and test set 
stratified_split('test_tsv.tsv', 'train_tsv', 'test_tsv', ratio=0.8)
#generate training vcf 
filter_vcf_by_samples('train.tsv', 'data\utils\ch_22_filtered.vcf', 'train.vcf')
#Generate test vcf 
filter_vcf_by_samples('test.tsv','data\utils\ch_22_filtered.vcf', 'test.vcf')



import h5py
import numpy as np
# import sys
# import glob
import allel
import pandas as pd


def vcf_to_npy(vcf_file, output_file):
    # Read VCF file
    vcf_data = allel.read_vcf(vcf_file)
    
    # Extract genotype data and transpose
    snps = vcf_data['calldata/GT'].transpose(1, 2, 0)  # Shape: (n_samples, n_alleles, n_variants)
    
    # Reshape the SNP data
    n_seq, n_chann, n_snps = snps.shape
    snps = snps.reshape(n_seq * n_chann, n_snps)  # Shape: (n_samples * n_alleles, n_variants)
    
    # Extract sample names
    samples = vcf_data['samples']
    
    # Save SNP data and sample names to .npy files
    np.save(f"{output_file}_snps.npy", snps)
    np.save(f"{output_file}_samples.npy", samples)

import numpy as np
import pandas as pd

def load_ancestry_labels(tsv_file, sample_names, output_file):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['sample_name', 'ancestry_label'])
    
    # Ensure that the sample names are in the same order as those in snps.npy
    # Match the sample names and create the ancestry labels array
    ancestry_labels = df.set_index('sample_name').loc[sample_names, 'ancestry_label'].values
    
    # Get unique ancestry labels and create a mapping to integer labels
    unique_labels = np.unique(ancestry_labels)
    label_mapping = {label: idx for idx, label in enumerate(unique_labels)}
    
    # Print the label mapping
    print("Label Mapping:")
    for label, idx in label_mapping.items():
        print(f"{label}: {idx}")
    
    # Replace ancestry labels with integer labels using the mapping
    ancestry_labels_int = np.array([label_mapping[label] for label in ancestry_labels])
    
    # Save the integer ancestry labels as a .npy file
    np.save(f"{output_file}_ancestry_labels.npy", ancestry_labels_int)
    
    return ancestry_labels_int


def generate_dataset(genotype_npy, labels_npy, output_file):
    # Load the combined genotype and label data from .npy files
    snps_data = np.load(genotype_npy)
    labels_data = np.load(labels_npy)

    # Create the output HDF5 file
    with h5py.File(output_file, "w") as outfile:
        # Create datasets for SNP data and labels
        outfile.create_dataset("vcf", data=snps_data)
        outfile.create_dataset("labels", data=labels_data)

    print(f"Dataset saved to {output_file}")


def vcf_to_npy(vcf_file, output_file):
    #Read VCF file
    vcf_data = allel.read_vcf(vcf_file)
    
    # Extract genotype data and transpose
    snps = vcf_data['calldata/GT'].transpose(1, 2, 0)  # Shape: (n_samples, n_alleles, n_variants)
    
    # Reshape the SNP data
    print("Reshaping the snp array")
    n_seq, n_chann, n_snps = snps.shape
    snps = snps.reshape(n_seq * n_chann, n_snps)  # Shape: (n_samples * n_alleles, n_variants)
    
    # Extract sample names
    print("Extracting samples")
    samples = vcf_data['samples']
    
    np.save(f"{output_file}_snps.npy", snps)
    np.save(f"{output_file}_samples.npy", samples)


def load_ancestry_labels(tsv_file, sample_names, output_file):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['sample_name', 'ancestry_label'])
    
    # Ensure that the sample names are in the same order as those in snps.npy
    # Match the sample names and create the ancestry labels array
    ancestry_labels = df.set_index('sample_name').loc[sample_names, 'ancestry_label'].values
    
    # Convert the ancestry labels into a 2D array (shape: [n_samples, 1])
    ancestry_labels_2d = ancestry_labels.reshape(-1, 1)  # This converts it to a 2D array
    
    # Save the ancestry labels as a .npy file
    np.save(f"{output_file}_ancestry_labels.npy", ancestry_labels_2d)
    
    return ancestry_labels_2d


def generate_dataset(genotype_npy, labels_npy, output_file):
    # Load the combined genotype and label data from .npy files
    snps_data = np.load(genotype_npy,allow_pickle=True)
    labels_data = np.load(labels_npy,allow_pickle=True)

    # Create the output HDF5 file
    with h5py.File(output_file, "w") as outfile:
        # Create datasets for SNP data and labels
        outfile.create_dataset("vcf", data=snps_data)
        outfile.create_dataset("labels", data=labels_data)

    print(f"Dataset saved to {output_file}")

# Example usage
if __name__ == '__main__':
    # Paths for input and output files
    vcf_file = r"E:\GeMorph\Ancestry\SALAI-Net\data\dataset_test.vcf"
    tsv_file = r"E:\GeMorph\Ancestry\SALAI-Net\data\dataset_test.tsv"
    output_file = r"E:\GeMorph\Ancestry\SALAI-Net\data\final_data\test"
    
    # Convert VCF to SNP and sample .npy files
    # vcf_to_npy(vcf_file, output_file)
    
    # Load ancestry labels and convert them to integer labels
    # snps = np.load(f"{output_file}_snps.npy")
    # sample_names = np.load(f"{output_file}_samples.npy",allow_pickle=True)
    # ancestry_labels = load_ancestry_labels(tsv_file, sample_names, output_file)


    # snps = np.load(r'E:\GeMorph\Ancestry\SALAI-Net\data\final_data\test_ancestry_labels.npy',allow_pickle=True)
    # print(snps)

    
    # Generate the HDF5 dataset from the .npy files
    # generate_dataset(f"{output_file}_snps.npy", f"{output_file}_ancestry_labels.npy", f"{output_file}_vcf_and_labels.h5")




# if __name__ == '__main__':

    # folder_name = sys.argv[1]

    # outfile = h5py.File(folder_name + "/vcf_and_labels.h5", "w")
    # # outfile = h5py.File("kk", "w")


    # gen_folder = glob.glob(folder_name + "/gen_*")
    # gen_folder.sort()

    # all_data = []
    # for gen in gen_folder:
    #     print(gen)
    #     all_data.append(np.load(gen + "/mat_vcf_2d.npy"))

    # all_data = np.concatenate(all_data)
    # outfile.create_dataset("vcf", data=all_data)

    # all_data = []
    # for gen in gen_folder:
    #     print(gen)
    #     all_data.append(np.load(gen + "/mat_map.npy"))

    # all_data = np.concatenate(all_data)
    # outfile.create_dataset("labels", data=all_data)
    # del all_data

    # outfile.close()


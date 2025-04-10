#!/usr/bin/env python

'''A function to get the sample data needed for the pipeline.'''

import pandas as pd
import glob
import os
import re
import sys

def get_sample_data(csv_file, input_dir):
    """
    A function to get the sample data needed for the pipeline.
    Supports multiple input formats:
    1. Standard Illumina paired-end fastq files with sample directories
    2. Single fastq inputs
    3. Paired-end fastqs with standard Illumina naming conventions
    """
    df = pd.read_csv(csv_file)

    # Assuming the sample names are in a column named 'sample'
    sample_names = df['sample'].astype(str).tolist()

    # Dictionary to hold the fastq file paths for each sample
    sample_fastq_files = {}
    all_samples_in_dir = set()
    
    # Check if samples are in subdirectories or directly in input_dir
    if any(os.path.isdir(os.path.join(input_dir, sample.split('_')[0])) for sample in sample_names):
        # Samples are in subdirectories
        all_samples_in_dir = set([d for d in os.listdir(input_dir) 
                                if os.path.isdir(os.path.join(input_dir, d))])
        use_subdirs = True
    else:
        # Samples might be directly in the input directory as fastq files
        all_samples_in_dir = set([s.split('_')[0] for s in os.listdir(input_dir) 
                                if s.endswith('.fastq.gz') or s.endswith('.fq.gz') or s.endswith('.fastq') or s.endswith('.fq')])
        use_subdirs = False

    missing_files = []
    not_found_samples = []

    for sample in sample_names:
        base_sample = sample.split('_')[0]  # Extract base sample name (without S*)
        sample_fastq_files[sample] = []
        
        if use_subdirs and base_sample not in all_samples_in_dir:
            not_found_samples.append(sample)
            continue
            
        for lane in range(1, 2):  # Loop through lanes 1 to 4
            if use_subdirs:
                # Standard format with sample directories
                patterns = [
                    # Standard Illumina format with .gz
                    (os.path.join(input_dir, base_sample, "{}_S*_L00{}_R1_*.fastq.gz".format(base_sample, lane)),
                     os.path.join(input_dir, base_sample, "{}_S*_L00{}_R2_*.fastq.gz".format(base_sample, lane)),
                     r'_S(\d+)_L00(\d)_R([12])_'),
                    # Standard Illumina format without .gz
                    (os.path.join(input_dir, base_sample, "{}_S*_L00{}_R1_*.fastq".format(base_sample, lane)),
                     os.path.join(input_dir, base_sample, "{}_S*_L00{}_R2_*.fastq".format(base_sample, lane)),
                     r'_S(\d+)_L00(\d)_R([12])_'),
                    # Simple paired format with .gz
                    (os.path.join(input_dir, base_sample, "{}*_R1*.fastq.gz".format(base_sample)),
                     os.path.join(input_dir, base_sample, "{}*_R2*.fastq.gz".format(base_sample)),
                     r'_R([12])'),
                    # Simple paired format without .gz
                    (os.path.join(input_dir, base_sample, "{}*_R1*.fastq".format(base_sample)),
                     os.path.join(input_dir, base_sample, "{}*_R2*.fastq".format(base_sample)),
                     r'_R([12])'),
                    # Single-end format with .gz
                    (os.path.join(input_dir, base_sample, "{}*.fastq.gz".format(base_sample)), 
                     None, None),
                    # Single-end format without .gz
                    (os.path.join(input_dir, base_sample, "{}*.fastq".format(base_sample)), 
                     None, None)
                ]
            else:
                # Files directly in input directory
                patterns = [
                    # Standard Illumina format with .gz
                    (os.path.join(input_dir, "{}_S*_L00{}_R1_*.fastq.gz".format(base_sample, lane)),
                     os.path.join(input_dir, "{}_S*_L00{}_R2_*.fastq.gz".format(base_sample, lane)),
                     r'_S(\d+)_L00(\d)_R([12])_'),
                    # Standard Illumina format without .gz
                    (os.path.join(input_dir, "{}_S*_L00{}_R1_*.fastq".format(base_sample, lane)),
                     os.path.join(input_dir, "{}_S*_L00{}_R2_*.fastq".format(base_sample, lane)),
                     r'_S(\d+)_L00(\d)_R([12])_'),
                    # Simple paired format with .gz
                    (os.path.join(input_dir, "{}*_R1*.fastq.gz".format(base_sample)),
                     os.path.join(input_dir, "{}*_R2*.fastq.gz".format(base_sample)),
                     r'_R([12])'),
                    # Simple paired format without .gz
                    (os.path.join(input_dir, "{}*_R1*.fastq".format(base_sample)),
                     os.path.join(input_dir, "{}*_R2*.fastq".format(base_sample)),
                     r'_R([12])'),
                    # Single-end format with .gz
                    (os.path.join(input_dir, "{}*.fastq.gz".format(base_sample)), 
                     None, None),
                    # Single-end format without .gz
                    (os.path.join(input_dir, "{}*.fastq".format(base_sample)), 
                     None, None)
                ]
            
            files_found = False
            
            for r1_pattern, r2_pattern, regex_pattern in patterns:
                r1_files = glob.glob(r1_pattern)
                r2_files = glob.glob(r2_pattern) if r2_pattern else []
                
                if r1_files:
                    files_found = True
                    
                    if regex_pattern and r'R([12])' in regex_pattern:
                        # Simple paired format with R1/R2
                        for file in r1_files:
                            sample_fastq_files[sample].append(('1', 'L00{}'.format(lane), 'R1', file))
                        for file in r2_files:
                            sample_fastq_files[sample].append(('1', 'L00{}'.format(lane), 'R2', file))
                    elif regex_pattern:
                        # Standard Illumina format - extract actual S number from filename
                        for file in r1_files:
                            match = re.search(regex_pattern, file)
                            if match:
                                # Use the actual sample number from the file instead of forcing S1
                                sample_number = match.group(1)
                                matched_lane = match.group(2)
                                sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R1', file))
                        
                        for file in r2_files:
                            match = re.search(regex_pattern, file)
                            if match:
                                # Use the actual sample number from the file instead of forcing S1
                                sample_number = match.group(1)
                                matched_lane = match.group(2)
                                sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R2', file))
                    else:
                        # Single-end format
                        for file in r1_files:
                            sample_fastq_files[sample].append(('1', 'L00{}'.format(lane), 'R1', file))
                    
                    # Break out of the patterns loop if files are found for this pattern
                    break
            
            if not files_found and lane == 1:
                missing_files.append("No fastq files found for sample {}".format(sample))
                break  # Skip checking other lanes if no files found for lane 1

    # Remove samples with no files
    empty_samples = [sample for sample, files in sample_fastq_files.items() if not files]
    for sample in empty_samples:
        del sample_fastq_files[sample]
        if sample not in not_found_samples:
            not_found_samples.append(sample)

    # Check for missing files
    if missing_files:
        for missing_file in missing_files:
            print(missing_file)
            
        # give option to continue without the missing files
        print("Do you want to continue without the missing files? (y/n)")
        response = input()
        if response.lower() == 'n':
            print("Pipeline terminated.")
            sys.exit(1)
        elif response.lower() == 'y':
            print("Continuing without the missing files.")
        else:
            raise ValueError("Invalid response. Please enter 'y' or 'n'.")
        
    # Print the number of found samples and samples not mentioned in the input list but present in the input directory
    found_samples = len(sample_fastq_files)
    extra_samples = all_samples_in_dir - set([s.split('_')[0] for s in sample_names])
    print(f"Found {found_samples} samples.")
    print(f"{len(extra_samples)} samples are not mentioned in the input list but present in the input directory: {', '.join(extra_samples)}")

    # Check for samples in the input list but not found in the input directory
    if not_found_samples:
        print(f"{len(not_found_samples)} samples are mentioned in the input list but not found in the input directory: {', '.join(not_found_samples)}")
        print("Do you want to start the analysis anyway? (y/n)")
        response = input()
        if response.lower() == 'n':
            print("Pipeline terminated.")
            sys.exit(1)
        elif response.lower() == 'y':
            print("Continuing without the missing samples.")
        else:
            raise ValueError("Invalid response. Please enter 'y' or 'n'.")

    # Merge the metadata with the fastq file information
    metadata_columns = df.columns.tolist()
    metadata_columns.remove('sample')

    output_rows = []

    for index, row in df.iterrows():
        sample = str(row['sample'])
        base_sample = sample.split('_')[0]  # Extract base sample name without S number
        metadata = [row[col] for col in metadata_columns]
        
        # Check if the sample exists in sample_fastq_files before proceeding
        if sample in sample_fastq_files:
            for sample_number, lane, read, file in sample_fastq_files[sample]:
                # Use the actual sample number from the file
                merged_sample = "{}_S{}".format(base_sample, sample_number)
                output_rows.append([merged_sample, lane, read, file] + metadata)

    # Save the results to a CSV file
    output_file = 'sample_data.csv'
    with open(output_file, 'w') as f:
        header = ['sample', 'lane', 'read', 'file'] + metadata_columns
        f.write(','.join(header) + '\n')
        
        for row in output_rows:
            f.write(','.join(map(str, row)) + '\n')
    
    # Return the CSV file as a dataframe
    return pd.read_csv(output_file)
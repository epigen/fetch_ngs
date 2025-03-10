#!/usr/bin/env python

#### libraries
import csv
import os, re, glob
import pandas as pd

#### configurations

# input
metadata_files = snakemake.input.metadata

# output
final_metadata_path = snakemake.output.metadata

# params
result_path = snakemake.params.result_path
output_fmt = snakemake.params.output_fmt      # "fastq" or "bam"
metadata_only = snakemake.params.metadata_only # 0 or 1

# set pu variables
metadata_tables = []
fieldnames = ["accession"]
accessions = []

# Process each input metadata file
for meta_file in metadata_files:
    # Derive the accession ID from the file name
    acc = os.path.basename(meta_file).replace(".metadata.csv", "")
    accessions.append(acc)
    # Read a snippet to detect delimiter
    with open(meta_file) as mf:
        sample_text = mf.read(1024)
        dialect = csv.Sniffer().sniff(sample_text.replace('\r\n', '\n'))
        delimiter = dialect.delimiter
    # Read the metadata CSV
    with open(meta_file, newline='') as mf:
        reader = csv.DictReader(mf, delimiter=delimiter)
        # Update the union of fieldnames
        for fn in reader.fieldnames:
            if fn not in fieldnames:
                fieldnames.append(fn)
        # For each row, add the accession information and collect it
        for row in reader:
            row["accession"] = acc
            metadata_tables.append(row)

# Create a combined metadata pandas df
metadata_df = pd.DataFrame(metadata_tables)
# replace empty strings with NA
metadata_df.replace(r'^\s*$', pd.NA, regex=True, inplace=True)
# drop empty columns
metadata_df.dropna(axis=1, how='all', inplace=True)

# if files have been downloaded add paths as second/third columns to df
if metadata_only==0:
    # Add file path columns based on output format
    if output_fmt == "fastq":
        suffix = "fastq.gz"
        pattern = re.compile(r'^(?P<sample>.+?)(?:_([12]))?\.fastq\.gz$')
        if "fastq_1" not in metadata_df.columns: metadata_df["fastq_1"] = pd.NA
        if "fastq_2" not in metadata_df.columns: metadata_df["fastq_2"] = pd.NA
    else:
        suffix = "bam"
        pattern = re.compile(r'^(?P<sample>.+?)\.bam$')
        if "bam" not in metadata_df.columns: metadata_df["bam"] = pd.NA
            
    # Find all FASTQ or BAM files across accessions and add them to the metadata
    for acc in accessions:
        acc_dir = os.path.join(result_path, acc)
        seq_files = glob.glob(os.path.join(acc_dir, f"*.{suffix}"))
        for f in seq_files:
            m = pattern.match(os.path.basename(f))
            if not m:
                continue
            
            sample = m.group("sample")
            pair = m.group(2) if output_fmt == "fastq" else None
            rows = metadata_df[metadata_df['accession'] == acc]
            matches = []
            
            for idx, row in rows.iterrows():
                if any(pd.notna(val) and sample in str(val) for val in row):
                    matches.append(idx)
            
            if len(matches) == 0:
                print(f"Error: No rows found for sample {sample} in accession {acc}. File: {f}")
                continue
            if len(matches) > 1:
                print(f"Error: Multiple rows found for sample {sample} in accession {acc}. File: {f}")
                continue
            row_index = matches[0]
            
            if output_fmt == "fastq":
                col = "fastq_1" if (pair is None or pair == "1") else "fastq_2"
                if pd.notna(metadata_df.at[row_index, col]):
                    print(f"Error: Duplicate assignment for sample {sample} in accession {acc} for column {col}. File: {f}")
                else:
                    metadata_df.at[row_index, col] = f
            else:
                col = "bam"
                if pd.notna(metadata_df.at[row_index, col]):
                    print(f"Error: Duplicate assignment for sample {sample} in accession {acc}. File: {f}")
                else:
                    metadata_df.at[row_index, col] = f

    # move file path columns to front
    if output_fmt == "fastq":
        file_cols = [col for col in ["fastq_1", "fastq_2"] if col in metadata_df.columns]
    else:
        file_cols = [col for col in ["bam"] if col in metadata_df.columns]
    other_cols = [col for col in metadata_df.columns if col not in file_cols]
    metadata_df = metadata_df[file_cols + other_cols]

# move accession to the first position
metadata_df = metadata_df[['accession'] + [col for col in metadata_df.columns if col != 'accession']]
# save final metadata df
metadata_df.to_csv(final_metadata_path, index=False)

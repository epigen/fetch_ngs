#!/usr/bin/env python
"""
merge_annotations.py

This script reads the metadata CSV files generated for each accession (by iSeq) and merges
them into one combined annotation CSV file for downstream analyses.

Steps:
1. For each metadata file (ending with .metadata.csv) in snakemake.input, detect the CSV dialect,
   and load its content.
2. Append the accession id as a column ("Accession") to every row.
3. Build a union of all column names encountered.
4. For each entry, add file path columns based on the desired output format (fastq or bam).
5. Write the merged annotation to the output file.
"""

import csv
import os

# Global variables provided by Snakemake
result_path = snakemake.params.result_path  # Base directory for results
output_fmt = snakemake.params.output_fmt      # "fastq" or "bam"
# accession_ids = snakemake.params.accession_ids  # List of accession IDs

metadata_tables = []
fieldnames = ["Accession"]

# Process each input metadata file
for meta_file in snakemake.input:
    # Only process metadata files, skip any marker files
    if meta_file.endswith(".metadata.csv"):
        # Derive the accession ID from the file name
        acc = os.path.basename(meta_file).replace(".metadata.csv", "")
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
                row["Accession"] = acc
                metadata_tables.append(row)

# Add file path columns based on output format
if output_fmt == "fastq":
    if "fastq_1" not in fieldnames:
        fieldnames.extend(["fastq_1", "fastq_2"])
else:
    if "bam" not in fieldnames:
        fieldnames.append("bam")

# For each metadata row, determine file paths for FASTQ or BAM files
for row in metadata_tables:
    acc = row["Accession"]
    # Construct the accession folder path
    acc_dir = os.path.join(result_path, acc)
    if output_fmt == "fastq":
        # Use Experiment or Run field if available; default to accession id
        exp_id = row.get("Experiment") or row.get("Run") or acc
        # Determine library layout and file naming conventions for paired or single-end data
        layout = (row.get("LibraryLayout") or row.get("Layout") or "").upper()
        if layout == "PAIRED" or os.path.exists(os.path.join(acc_dir, f"{exp_id}_2.fastq.gz")):
            row["fastq_1"] = os.path.join(acc_dir, f"{exp_id}_1.fastq.gz")
            row["fastq_2"] = os.path.join(acc_dir, f"{exp_id}_2.fastq.gz")
        else:
            row["fastq_1"] = os.path.join(acc_dir, f"{exp_id}.fastq.gz")
            row["fastq_2"] = ""
    else:
        # For BAM output, assume the BAM file is named with the experiment id or accession id
        exp_id = row.get("Experiment") or row.get("Run") or acc
        row["bam"] = os.path.join(acc_dir, f"{exp_id}.bam")

# Write out the combined annotation CSV file
with open(snakemake.output.annotation, 'w', newline='') as out_csv:
    writer = csv.DictWriter(out_csv, fieldnames=fieldnames)
    writer.writeheader()
    for row in metadata_tables:
        writer.writerow(row)

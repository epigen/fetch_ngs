#!/usr/bin/env python
"""
convert_fastq_to_bam.py

This script scans an accession folder for FASTQ files (which may have arbitrary names),
groups them into samples based on a naming pattern (e.g. detecting _1 and _2 for paired-end),
and then converts each group to an unmapped BAM file using Picard’s FastqToSam.

Steps:
1. List all FASTQ (.fastq.gz) files in the provided accession directory.
2. Group files by a common sample prefix:
   - If a file name matches the pattern (.*)_[12].fastq.*.gz, the captured group is used as the sample name.
   - Otherwise, the file name (without extension) is used as the sample name.
3. For each sample group:
   - If there are two files, treat them as paired-end (F1 and F2).
   - If there is only one file, treat it as single-end.
   - Invoke Picard FastqToSam with appropriate arguments to produce {sample}.bam.
4. Remove the FASTQ files after successful conversion.
"""

import os
import glob
import re
import subprocess
import argparse

def group_fastq_files(acc_dir):
    """Group FASTQ files in acc_dir by sample name."""
    fastq_files = glob.glob(os.path.join(acc_dir, "*.fastq.gz"))
    groups = {}
    for fq in fastq_files:
        fname = os.path.basename(fq)
        # Check for paired-end pattern like sample_1.fastq.gz or sample_R1.fastq.gz
        m = re.match(r"(.*?)(?:_[12]|_R[12])\.fastq.*\.gz", fname, re.IGNORECASE)
        if m:
            sample = m.group(1)
        else:
            # Fallback: remove the extension and use the whole name as sample
            sample = os.path.splitext(os.path.splitext(fname)[0])[0]
        groups.setdefault(sample, []).append(fq)
    return groups

def convert_group_to_bam(sample, files, acc_dir):
    """Convert a group of FASTQ files to BAM using Picard FastqToSam."""
    bam_output = os.path.join(acc_dir, f"{sample}.bam")
    cmd = ["picard", "FastqToSam", f"OUTPUT={bam_output}",
           f"SAMPLE_NAME={sample}", f"READ_GROUP_NAME={sample}"]
    # If two files are detected, assume paired-end
    if len(files) == 2:
        # Sort files so that F1 comes before F2 if not already ordered
        files.sort()
        cmd.extend([f"F1={files[0]}", f"F2={files[1]}"])
    elif len(files) == 1:
        cmd.append(f"F1={files[0]}")
    else:
        raise ValueError(f"Sample '{sample}' has {len(files)} FASTQ files; expected 1 (single-end) or 2 (paired-end).")
    print(f"Converting sample '{sample}' with files: {files}")
    subprocess.check_call(cmd)
    return bam_output

def remove_fastq_files(files):
    """Remove the given FASTQ files."""
    for f in files:
        os.remove(f)

def main():
    parser = argparse.ArgumentParser(description="Convert FASTQ files to unmapped BAM using Picard.")
    parser.add_argument("--acc_dir", required=True, help="Path to accession directory containing FASTQ files")
    args = parser.parse_args()
    acc_dir = args.acc_dir

    groups = group_fastq_files(acc_dir)
    if not groups:
        print(f"No FASTQ files found in {acc_dir}.")
        return

    for sample, files in groups.items():
        try:
            convert_group_to_bam(sample, files, acc_dir)
            remove_fastq_files(files)
        except Exception as e:
            print(f"Error processing sample '{sample}': {e}")
            raise

if __name__ == "__main__":
    main()

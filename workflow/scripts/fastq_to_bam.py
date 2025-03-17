#!/usr/bin/env python

#### libraries
import os, re, glob, subprocess

#### configurations

# output
marker_file = snakemake.output.marker

# params
acc_dir = snakemake.params.acc_dir

# identify all FASTQ files
fastq_files = glob.glob(os.path.join(acc_dir, "*.fastq.gz"))
samples = {}
pattern = re.compile(r"(.*?)(?:_([12]))?\.fastq\.gz$")

# identify sample names and type (single/paired-end) across files
for f in fastq_files:
    m = pattern.match(os.path.basename(f))
    if not m:
        continue
    sample, pair = m.group(1), m.group(2)
    samples.setdefault(sample, {})
    if pair:
        samples[sample][f"r{pair}"] = f
    else:
        samples[sample]["se"] = f

bam_files = []  # list to track processed .bam files

# build and execute Picard FastqToSam command depending on type
for sample, fdict in samples.items():
    outbam = os.path.join(acc_dir, f"{sample}.bam")
    
    if os.path.exists(outbam):  # Skip if outbam already exists
        bam_files.append(outbam)
        continue

    if "r1" in fdict and "r2" in fdict:
        cmd = (
            f"picard FastqToSam FASTQ={fdict['r1']} FASTQ2={fdict['r2']} "
            f"OUTPUT={outbam} SAMPLE_NAME={sample}"
        )
    elif "se" in fdict:
        cmd = (
            f"picard FastqToSam FASTQ={fdict['se']} "
            f"OUTPUT={outbam} SAMPLE_NAME={sample}"
        )
    else:
        continue

    subprocess.run(cmd, shell=True, check=True)
    bam_files.append(outbam)  # log the newly created .bam

# log created bam files as success
if len(bam_files)>0:
    with open(marker_file, "w") as f:
        for bam in bam_files:
            f.write(f"{bam}\n")
else:
    print(f"Error: No uBAM samples created or found in {acc_dir}.")

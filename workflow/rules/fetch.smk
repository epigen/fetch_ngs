
# download metadata and fastq.gz files (optional) for each accession
rule iseq_download:
    output:
        acc_dir = directory(os.path.join(result_path, "{accession}")),
        metadata = os.path.join(result_path, "{accession}", "{accession}.metadata.csv"),
        success_log = os.path.join(result_path, "{accession}", "success.log") if config["metadata_only"]==0 else [],
    params:
        accession = lambda wc: "{}".format(wc.accession),
        metadata_only = "-m" if config["metadata_only"]==1 else "",
    resources:
        mem_mb=config.get("mem", "16000"),
        parallel_downloads=1, # https://stackoverflow.com/questions/51977436/restrict-number-of-jobs-by-a-rule-in-snakemake
    threads: config.get("threads", 4)
    conda:
        "../envs/iseq.yaml"
    shell:
        """
        iseq -i {params.accession} -o {output.acc_dir} -g -t {threads} -p {threads} {params.metadata_only}
        # if only a TSV metadata file is generated, convert it to CSV using sed (replace tabs with commas)
        if [ -f {output.acc_dir}/{params.accession}.metadata.tsv ]; then \
            awk 'BEGIN {{ FS="\\t"; OFS="," }} {{
              rebuilt=0
              for(i=1;i<=NF;++i) {{
                if ($i ~ /,/ && $i !~ /^".*"$/) {{
                  gsub("\\"", "\\"\\"", $i)
                  $i = "\\"" $i "\\""
                  rebuilt=1
                }}
              }}
              if (!rebuilt) {{ $1=$1 }}
              print
            }}' {output.acc_dir}/{params.accession}.metadata.tsv > {output.acc_dir}/{params.accession}.metadata.csv; \
        fi
        """

# optional: convert fastq.gz to uBAM using Picard
rule fastq_to_bam:
    input:
        metadata = os.path.join(result_path, "{accession}", "{accession}.metadata.csv")
    output:
        marker = os.path.join(result_path,".fastq_to_bam","{accession}.done")
    params:
        acc_dir = lambda wc: os.path.join(result_path, "{}".format(wc.accession)),    
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/picard.yaml"
    script:
        "../scripts/fastq_to_bam.py"

# phantom/shadow rule to enable downstream processing
# requires to know that the file will exist in that exact location, otherwise MissingInputException error
rule fetch_file:
    input:
        metadata = os.path.join(result_path, "{accession}", "{accession}.metadata.csv"),
        bam_confirmation = os.path.join(result_path, ".fastq_to_bam","{accession}.done") if output_fmt=="bam" else [],
    output:
        seqfile = update(os.path.join(result_path,"{accession}",r"{sample}.{suffix,(bam|fastq\.gz)}")),
    params:
        acc_dir = lambda wc: os.path.join(result_path, "{}".format(wc.accession)),    
    resources:
        mem_mb="1000",
    shell:
        """
        # only if the file already exists
        if [ -f {output.seqfile} ]; then \
            touch {output.seqfile}; \
        fi
        """
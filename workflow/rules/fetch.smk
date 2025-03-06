
rule iseq_download:
    output:
        acc_dir = directory(os.path.join(result_path, "{accession}")),
        metadata = os.path.join(result_path, "{accession}", "{accession}.metadata.csv")
    params:
        accession = lambda wc: "{}".format(wc.accession),
    threads: 8
    conda:
        "../envs/iseq.yaml"
    shell:
        """
        iseq -i {params.accession} -o {output.acc_dir} -g -t {threads}
        # if a TSV metadata file is generated, convert it to CSV using sed (replace tabs with commas)
        if [ -f {output.acc_dir}/{params.accession}.metadata.tsv ]; then \
            sed 's/\\t/,/g' {output.acc_dir}/{params.accession}.metadata.tsv > {output.acc_dir}/{params.accession}.metadata.csv; \
        fi
        """

# rule fastq_to_bam:
#     input:
#         acc_dir = os.path.join(result_path, "{accession}")
#     output:
#         marker = os.path.join(result_path, "{accession}", ".fastq_to_bam.done")
#     conda:
#         "../envs/picard.yaml"
#     shell:
#         """
#         cd {input.acc_dir} && \
#         if [ -f "{wildcards.accession}_1.fastq.gz" ]; then \
#           if [ -f "{wildcards.accession}_2.fastq.gz" ]; then \
#             picard FastqToSam F1={wildcards.accession}_1.fastq.gz F2={wildcards.accession}_2.fastq.gz OUTPUT={wildcards.accession}.bam SAMPLE_NAME={wildcards.accession} READ_GROUP_NAME={wildcards.accession}; \
#           else \
#             picard FastqToSam F1={wildcards.accession}_1.fastq.gz OUTPUT={wildcards.accession}.bam SAMPLE_NAME={wildcards.accession} READ_GROUP_NAME={wildcards.accession}; \
#           fi; \
#         fi; \
#         rm -f {wildcards.accession}_1.fastq.gz {wildcards.accession}_2.fastq.gz; \
#         touch {output.marker}
#         """
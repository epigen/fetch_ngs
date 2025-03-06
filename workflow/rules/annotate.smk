
rule merge_annotation:
    input:
        expand(os.path.join(result_path, "{accession}", "{accession}.metadata.csv"), accession = accession_ids),
        expand(os.path.join(result_path, "{accession}", ".fastq_to_bam.done"), accession = accession_ids) if output_fmt=="bam" else [],
    output:
        annotation = os.path.join(result_path, "annotation.csv")
    params:
        result_path = result_path,
        output_fmt = output_fmt,
        accession_ids = accession_ids
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_annotations.py"

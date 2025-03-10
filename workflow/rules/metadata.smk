
# merge accession-wise metadata into one comprehensive metadata file 
rule merge_metadata:
    input:
        metadata = expand(os.path.join(result_path, "{accession}", "{accession}.metadata.csv"), accession = accession_ids),
        bam_confirmation = expand(os.path.join(result_path, ".fastq_to_bam","{accession}.done"), accession = accession_ids) if output_fmt=="bam" and config["metadata_only"]==0 else [],
    output:
        metadata = report(os.path.join(result_path, "metadata.csv"),
                                caption="../report/metadata.rst",
                                category="{}_{}".format(config["project_name"], module_name),
                                subcategory="Metadata",
                                labels={
                                    "name": "Metadata",
                                    "type": "CSV",
                                }),
    params:
        result_path = result_path,
        output_fmt = output_fmt,
        metadata_only = config["metadata_only"],
    resources:
        mem_mb="4000",
    conda:
        "../envs/picard.yaml"
    script:
        "../scripts/merge_metadata.py"



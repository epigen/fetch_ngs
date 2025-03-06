def get_annotation_inputs(wc):
    files = []
    for acc in dataset_accessions[wc.dataset]:
        files.append(f"{result_path}/{wc.dataset}/fetch_ngs/{acc}/{acc}.metadata.csv")
        if output_fmt == "bam":
            files.append(f"{result_path}/{wc.dataset}/fetch_ngs/{acc}/.fastq_to_bam.done")
    return files
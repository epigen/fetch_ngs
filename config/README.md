You only need one configuration file to run the complete workflow. You can use the provided example as starting point. If in doubt read the comments in the configuration file, the documentation of the respective methods and/or try the default values.

**configuration (`config/config.yaml`):** Different for every project/dataset and configures the datasets to be fetched and how they should be processed. The fields are described within the file.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default/config.yaml`, which supersedes global Snakemake profiles.


## Example Configurations

### Metadata-only

```yaml
project_name: ExploratoryProject
result_path: results/

metadata_only: 1

accession_ids:
  - GSE122139
```

### Full download with BAM output

```yaml
threads: 16
mem: 32000

project_name: BAMProject
result_path: results/

metadata_only: 0

output_format: bam

accession_ids:
  - GSE122139
  - SRP123456
```

### Full download with FASTQ output

```yaml
threads: 16
mem: 32000

project_name: FastqProject
result_path: results/

metadata_only: 0

output_format: fastq

accession_ids:
  - ERS5684710
```

# ALS_RNAseq

Scripts for building the ALS bulk RNA-seq dataset: aligning and quantifying
FASTQs with nf-core/rnaseq via ORCA/Tower, with QC and expression
normalization stages planned next. Outputs are documented and versioned on
[Synapse](https://www.synapse.org/) for use by experimental biologists,
clinicians, and platform engineers.

See [`docs/methods.md`](docs/methods.md) for the full pipeline description
and reproducibility notes.

## Repository Structure

```
ALS_RNAseq/
├── docs/
│   └── methods.md                        # pipeline description, reproducibility notes
├── config/
│   ├── als_dataset_syn_ids.yml            # Synapse IDs per dataset (samplesheet, output folder)
│   └── pipeline_config.yml                # AWS profile, S3 bucket, per-dataset staging keys
├── scripts/
│   └── 01_alignment/                     # nf-core/rnaseq alignment + quantification
│       ├── orca_rnaseq_workflow.py        # orchestrates the ORCA/Tower pipeline
│       ├── split_merged_matrics_by_study.py # splits merged RSEM matrices per study, uploads to Synapse
│       ├── als_rnaseq_params.json        # nf-core/rnaseq science params
│       └── nextflow-rnaseq.config        # Nextflow/Tower infra config
│   # 02_qc/ and 03_normalization/ planned
├── samplesheets/
│   ├── dataset_targets.csv
│   ├── samplesheet.csv
│   └── samplesheet_with_dataset_names.csv
└── data/                                  # local working files (gitignored, not committed)
    ├── rsem_merged_tsv/                   # merged RSEM matrices downloaded from the pipeline run
    └── per_study/                         # per-study matrices produced by the split step
```

Scripts are grouped by pipeline stage (`01_alignment`, with `02_qc` and
`03_normalization` to follow) so the processing order is visible directly in
the file tree.

## Pipeline

1. **`scripts/01_alignment/orca_rnaseq_workflow.py`** — orchestrates the
   ORCA/Tower run in four steps (fetch samplesheet from Synapse, stage
   FASTQs, run nf-core/rnaseq, index results back to Synapse). Reads
   Synapse IDs and S3/AWS settings from `config/als_dataset_syn_ids.yml` and
   `config/pipeline_config.yml` -- add a cohort/disease by editing those
   files, not the script. See the script's docstring for setup and usage.
2. **`scripts/01_alignment/split_merged_matrics_by_study.py`** — once the run
   produces merged RSEM matrices, streams them into one file per study
   (`split`) and uploads the per-study files to Synapse (`upload`). See the
   script's docstring for usage, including `--check` and `--test` for
   verifying the sample-to-study mapping and split correctness without a
   full run.

Planned next: `scripts/02_qc/` (post-alignment QC) and
`scripts/03_normalization/` (expression normalization).

## Samplesheets

### `dataset_targets.csv`

Targets file listing the ALS bulk RNA-seq datasets selected from the Synapse
DatasetCollection, along with their Synapse dataset IDs and targeted FASTQ
subfolder IDs.

Columns: `dataset_name`, `dataset_synid`, `fastq_folder_synid`

### `samplesheet.csv`

nf-core/rnaseq samplesheet with FASTQ files referenced by Synapse ID.

Columns: `sample`, `single_end`, `fastq_1`, `fastq_2`, `strandedness`

### `samplesheet_with_dataset_names.csv`

Same as `samplesheet.csv` with an additional `dataset` column identifying the
source dataset for each sample. Used by `split_merged_matrics_by_study.py` to
map samples back to their study.

Columns: `sample`, `single_end`, `fastq_1`, `fastq_2`, `strandedness`, `dataset`

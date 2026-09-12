# Methods

Reference for collaborators (experimental biologists, clinicians, platform
engineers) and for methods sections describing how this ALS bulk RNA-seq
dataset was produced. Written to be reproducible by someone outside this
repo's original authors, given Synapse access and the credentials listed in
each stage's script.

## Pipeline overview

```
Synapse DatasetCollection (raw FASTQs)
        │
        ▼
┌───────────────────────┐
│ 01_alignment           │  nf-core/rnaseq (STAR + RSEM) via ORCA/Tower,
│                        │  GRCh38 / GENCODE v43. Produces merged RSEM
│                        │  matrices, split per study, uploaded to Synapse.
└───────────────────────┘
        │
        ▼
┌───────────────────────┐
│ 02_qc                  │  planned
└───────────────────────┘
        │
        ▼
┌───────────────────────┐
│ 03_normalization       │  planned
└───────────────────────┘
```

## Configuration

Synapse IDs and S3/AWS settings live in `config/als_dataset_syn_ids.yml` and
`config/pipeline_config.yml`, not in the scripts -- so the provenance record
(exactly which Synapse entity backs each input/output) is a reviewable config
diff rather than a code change, and new cohorts/diseases can be added without
editing Python.

- **`config/als_dataset_syn_ids.yml`** — one entry per dataset under `datasets`: the
  samplesheet's Synapse ID and filename, and the Synapse folder ID that
  `synindex` writes results back into.
- **`config/pipeline_config.yml`** — the AWS profile and S3 bucket used for
  staging, and each dataset's S3 staging key prefix.
- **`scripts/01_alignment/als_rnaseq_params.json`** — nf-core/rnaseq science
  parameters (aligner, skip flags, reference genome paths) shared across
  datasets; not dataset-specific, so it stays separate from `config/`.

To add a cohort/disease: append an entry to `config/als_dataset_syn_ids.yml`'s
`datasets` list and a matching key under `pipeline_config.yml`'s
`staging_keys`. `orca_rnaseq_workflow.py`'s `generate_datasets()` reads both
files directly -- no other code changes needed.

## 01_alignment

- **Pipeline:** [nf-core/rnaseq](https://nf-co.re/rnaseq), aligner `star_rsem`.
- **Reference:** GRCh38 primary assembly, GENCODE v43 annotation
  (see `scripts/01_alignment/als_rnaseq_params.json` for the exact
  `fasta`/`gtf` S3 paths and all other pipeline parameters).
- **Orchestration:** [ORCA](https://github.com/Sage-Bionetworks-Workflows/py-orca)
  against Seqera/Nextflow Tower, implementing
  [SKG-145](https://sagebionetworks.jira.com/browse/SKG-145). Four steps:
  fetch samplesheet -> stage FASTQs from Synapse -> run nf-core/rnaseq ->
  index results back to Synapse. See
  `scripts/01_alignment/orca_rnaseq_workflow.py` docstring for exact commands,
  required credentials, and environment variables.
- **Post-processing:** the pipeline's output is one merged matrix per
  quantification type (gene/transcript x counts/TPM) across all samples in
  the run. `scripts/01_alignment/split_merged_matrics_by_study.py` splits
  each merged matrix into one file per study (using
  `samplesheets/samplesheet_with_dataset_names.csv` as the sample-to-study
  map) and uploads the per-study files to Synapse, one folder per study.

## 02_qc (planned)

Not yet implemented.

## 03_normalization (planned)

Not yet implemented.

## Reproducibility notes

- Every stage's script docstring lists its own prerequisites (credentials,
  packages, environment variables) and exact CLI usage -- treat the
  docstring as the source of truth; this file gives the cross-stage picture.
- `scripts/01_alignment/split_merged_matrics_by_study.py --check` and
  `--test` verify the sample-to-study mapping and split correctness against
  real data without requiring a full pipeline run -- run these after editing
  `samplesheets/dataset_targets.csv` or the sample sheet.
- Large intermediate files (merged RSEM matrices, per-study splits) are
  staged locally under `data/` and are not committed to this repository
  (see `.gitignore`); the durable copies live on Synapse. See each stage's
  upload step for destination Synapse IDs.

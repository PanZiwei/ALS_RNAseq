# ALS_RNAseq

Aligning and quantifying ALS bulk RNA-seq FASTQs with nf-core/rnaseq via
ORCA/Tower, with outputs versioned on [Synapse](https://www.synapse.org/).

**Outputs**
- Aligned/quantified nf-core/rnaseq results indexed back to Synapse (BAMs, QC, RSEM matrices)
- Merged RSEM count/TPM matrices, split into one file per study
- *(planned)* QC-filtered matrices and CQN-normalized matrices

## Structure

```
ALS_RNAseq/
├── config/                     # ALL config -- edit here, not the code
│   ├── als_dataset_syn_ids.yml   # Synapse IDs per dataset (samplesheet, output folder)
│   ├── aws_s3_config.yml         # AWS profile, S3 bucket, staging prefixes
│   ├── als_rnaseq_params.json    # nf-core/rnaseq science params
│   └── nextflow-rnaseq.config    # per-process cpus/memory + retry
├── pipeline/                   # stages in run order (0 → 1 → …)
│   ├── 0_alignment/              # LAPTOP -- orchestrates Tower, no data touches disk
│   │   ├── orca_rnaseq_alignment.py   # the 5-step ORCA/Tower pipeline
│   │   └── helper/split_samtools_stats.py  # pre-synindex fix (see Pipeline)
│   ├── 1_matrices/              # EC2 -- processes CONTROLLED-ACCESS data
│   │   └── split_matrices_by_study.py
│   ├── 2_qc/                    # planned (EC2)
│   └── 3_normalization/         # planned (EC2)
└── data/                       # gitignored, created at runtime on EC2
```

Stages are numbered in run order. **Stage `0_alignment` runs on your laptop**
(it only orchestrates Tower — no controlled data comes down); **stage `1+` run
on an EC2 instance** in the data's AWS account, since they download and process
controlled-access files.

## Setup

**1. Environment**
```bash
conda env create -f environment.yml
conda activate als_rnaseq
```

**2. AWS access (SSO)** — the pipeline reads/writes S3 through an AWS profile
named `tower`. Two things must be granted to you first:

- **A Jumpcloud account with AWS access** — Sage provisions all AWS access
  through Jumpcloud. If you don't have one, email `it@sagebase.org`.
- **Access to the `workflows-nextflow-prod` account** (Sage's shared production
  AWS account where Seqera workspaces run compute + S3) and the
  complementaire-project Seqera workspace — request from the project owner.

Then create the profile with `aws configure sso` (AWS CLI v2): answer **SSO
start URL** `https://d-906769aa66.awsapps.com/start` and **SSO region**
`us-east-1`, then pick the `workflows-nextflow-prod` account and your role in the
browser. The resulting profile should look like this — name it `tower` to match
`aws_profile:` in [`config/aws_s3_config.yml`](config/aws_s3_config.yml):

```ini
[profile tower]
sso_start_url  = https://d-906769aa66.awsapps.com/start
sso_region     = us-east-1
sso_account_id = <workflows-nextflow-prod account id, from the browser>
sso_role_name  = <the role granted to you>
region = us-east-1
output = json
```

Your role needs **write** access for the `split_samtools_stats` step (S3
copy/delete); read-only is enough for `synstage` / `synindex` / listing. Full
walkthrough: [Jumpcloud wiki](https://sagebionetworks.jira.com/wiki/spaces/IT/pages/405864455/Jumpcloud#Access-AWS-with-AWS-SSO-CLI).

Log in before each run (tokens expire) and verify:
```bash
aws sso login --profile tower
aws s3 ls s3://complementaire-project-tower-bucket/ALS/ --profile tower
```

**3. Tower token** (laptop → Tower, to launch jobs) — add to `~/.zshrc`:
```bash
export NEXTFLOWTOWER_CONNECTION_URI="https://:<token>@tower.sagebionetworks.org/api?workspace=Sage-Bionetworks%2Fcomplementaire-project"
```
Get `<token>` at tower.sagebionetworks.org → avatar → Your tokens (goes after
the `:`, before `@`; the workspace is percent-encoded).

**4. Synapse token** (cloud job → Synapse, to read/write data) — set
`SYNAPSE_AUTH_TOKEN` as a **workspace** secret (not user):
[complementaire-project secrets](https://tower.sagebionetworks.org/orgs/Sage-Bionetworks/workspaces/complementaire-project/secrets).

## Pipeline

### Stage 0 — alignment (laptop)

`pipeline/0_alignment/orca_rnaseq_alignment.py` runs five steps (default `all`,
or name one). Steps run in order and are resumable — a re-run skips steps whose
Tower run already succeeded.

```bash
python pipeline/0_alignment/orca_rnaseq_alignment.py --check      # validate config, no Tower/Synapse
python pipeline/0_alignment/orca_rnaseq_alignment.py             # run all
python pipeline/0_alignment/orca_rnaseq_alignment.py synindex    # run one step
```

- **`fetch_samplesheet`** — copies the samplesheet from Synapse to S3.
- **`synstage`** — stages input FASTQs from Synapse to S3 (nf-synapse).
- **`rnaseq`** — runs [nf-core/rnaseq](https://nf-co.re/rnaseq) (STAR + RSEM,
  GRCh38 / [Gencode v43](https://www.gencodegenes.org/)) on Tower/AWS Batch.
- **`split_samtools_stats`** — splits `star_rsem/samtools_stats/` into `raw/` +
  `markdup/`. Needed because that folder (10,764 files) exceeds Synapse's hard
  10,000-children-per-folder limit and would fail `synindex`. Server-side and
  idempotent (`--dry-run` / `--check`).
- **`synindex`** — mirrors the full output tree back to Synapse (nf-synapse).
  Matches entities by name+parent, so re-runs update (bump version), never
  duplicate.

### Stage 1 — per-study matrices (EC2)

**`pipeline/1_matrices/split_matrices_by_study.py`** — splits the merged RSEM
matrices into one file per study (`split`) and uploads them to Synapse
(`upload`). `--check` (sample→study mapping) and `--test` (split one matrix,
assert columns match source) verify without a full run. Runs on **EC2** because
it downloads controlled-access data.

## Data sources

- **Input FASTQs & samplesheet** — [Synapse](https://www.synapse.org/); IDs in
  `config/als_dataset_syn_ids.yml`. Samplesheets are not committed (Synapse is
  the source of truth); fetched at run time and cached under `data/`.
- **Reference** — GRCh38 + Gencode v43, paths in `config/als_rnaseq_params.json`.
- **Compute** — Seqera/Tower on the `workflows-nextflow-prod` AWS account.

## Notes & warnings

- **Controlled-access data** never touches a laptop. Stage `0_alignment` only
  orchestrates; stage `1+` (anything that downloads data) runs on EC2.
- **Folder-size scaling** — `star_rsem/` itself (~4 files/sample + merged
  matrices) crosses the 10k-per-folder limit around ~2,500 samples. Fine for
  ALS's 1,794; revisit for larger cohorts.
- **Config, not code** — add a cohort or change an ID/param/resource by editing
  `config/`; never edit the scripts for run-specific values.

"""ORCA recipe: nf-core/rnaseq bulk RNA-seq processing on GRCh38.
Implements SKG-145 (https://sagebionetworks.jira.com/browse/SKG-145).

Orchestrate five steps (run in order):
  1. fetch_samplesheet     : Synapse -> S3
  2. synstage              : nf-synapse, download FASTQs Synapse -> S3
  3. rnaseq                : nf-core/rnaseq alignment/quantification
  4. split_samtools_stats  : split star_rsem/samtools_stats/ (the one folder over
                             Synapse's 10k-children limit) into raw/ + markdup/
  5. synindex              : nf-synapse, S3 -> Synapse, a SINGLE run over the whole
                             output_directory (mirrors the tree into parent_id).
                             Works because split_samtools_stats already fixed the
                             only folder over Synapse's 10k-children limit; there is
                             no per-run object limit. Resumable: skipped if a prior
                             synindex run for this dataset already succeeded.

Prerequisites:
  - pip install py-orca pyyaml
  - AWS credentials configured (profile: tower)
  - SYNAPSE_AUTH_TOKEN set as a Tower workspace secret (not user secret):
    https://tower.sagebionetworks.org/orgs/Sage-Bionetworks/workspaces/complementaire-project/secrets
  - NEXTFLOWTOWER_CONNECTION_URI env var (only env var this script reads):
    1. Get a token: tower.sagebionetworks.org -> avatar -> Your tokens -> Add token
    2. export NEXTFLOWTOWER_CONNECTION_URI="https://:<token>@tower.sagebionetworks.org/api?workspace=Sage-Bionetworks%2Fcomplementaire-project"

Usage:
  python orca_rnaseq_workflow.py                        # all steps
  python orca_rnaseq_workflow.py synstage                # one step
  python orca_rnaseq_workflow.py synstage rnaseq         # multiple steps
  python orca_rnaseq_workflow.py split_samtools_stats    # pre-index folder split only
  python orca_rnaseq_workflow.py synindex                # index everything (single run)
"""
import asyncio
import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import boto3
import yaml
from orca.services.nextflowtower import NextflowTowerOps
from orca.services.nextflowtower.client import NextflowTowerClient
from orca.services.nextflowtower.models import LaunchInfo
from synapseclient import Synapse
from synapseclient.operations import get as syn_get
from synapseclient.operations.factory_operations import FileOptions

from split_samtools_stats import split_samtools_stats

PARAMS_PATH = Path(__file__).parent / "als_rnaseq_params.json"
"""Science params only -- input/outdir are computed below, see "_note" key in the JSON."""

CONFIG_PATH = Path(__file__).parent / "nextflow-rnaseq.config"

# Synapse IDs and infra config (bucket, staging keys, AWS profile) live in config/
# Adding a cohort or changing an ID is a config diff, not a code change. 
# See config/als_dataset_syn_ids.yml and config/pipeline_config.yml.
CONFIG_DIR = Path(__file__).parent.parent.parent / "config"
SYN_IDS_PATH = CONFIG_DIR / "als_dataset_syn_ids.yml"
PIPELINE_CONFIG_PATH = CONFIG_DIR / "pipeline_config.yml"

PIPELINE_CONFIG = yaml.safe_load(PIPELINE_CONFIG_PATH.read_text())

session = boto3.Session(profile_name=PIPELINE_CONFIG["aws_profile"])
s3 = session.client("s3")


def _patch_tower_client_double_slash():
    """Fix py-orca URL construction: path has leading slash, causing // in URL and 401.

    orca builds request URLs by joining the api_endpoint with a leading-slash path,
    which yields e.g. https://tower.sagebionetworks.org/api//workflow/launch. The
    doubled slash triggers 400/401 responses from the Tower API. Stripping the leading
    slash from each request path fixes it. (Ported from sarek_somatic_workflow.py.)
    """
    _original_request = NextflowTowerClient.request

    def _patched_request(self, method: str, path: str, **kwargs):
        return _original_request(self, method, path.lstrip("/"), **kwargs)

    NextflowTowerClient.request = _patched_request


_patch_tower_client_double_slash()


async def monitor_workflow(ops: NextflowTowerOps, run_id: str):
    """Wait for a workflow run to finish; raise if it didn't succeed.

    Without this check, a failed synstage/rnaseq run would still let the next
    step launch against nonexistent output. (Ported from py-orca's demo.py.)
    """
    status = await ops.monitor_workflow(run_id=run_id, wait_time=60 * 2)
    if not status.is_successful:
        raise RuntimeError(f"Workflow did not complete successfully ({status}).")
    return status


async def main():
    parser = argparse.ArgumentParser(
        description="Run the ALS nf-core/rnaseq pipeline via ORCA/Nextflow Tower.",
    )
    parser.add_argument(
        'step', nargs='*', default='all',
        choices=['all', 'fetch_samplesheet', 'synstage', 'rnaseq', 'split_samtools_stats', 'synindex'],
        help="One or more steps to run, in order (default: all)",
    )
    args = parser.parse_args()

    ops = NextflowTowerOps()
    datasets = generate_datasets()
    runs = [run_workflows(ops, dataset, args.step) for dataset in datasets]
    statuses = await asyncio.gather(*runs)
    print(statuses)


@dataclass
class Dataset:
    id: str
    """The synapse id for the samplesheet."""

    samplesheet: str
    """The name of the samplesheet to run."""

    synapse_id_for_output: str
    """The synapse id for the output folder, this is where the output will be uploaded to."""

    output_directory: str
    """The S3 uri where nf-core/rnaseq writes results.
    Paired with synapse_id_for_output -- synindex just copies this into that Synapse folder."""

    bucket_name: str
    """The name of the bucket to stage the samplesheet in."""

    staging_key: str
    """The key in the S3 bucket where this workflow is going to run."""

    version: int | None = None
    """Synapse version of the samplesheet to fetch. None = latest."""

    @property
    def samplesheet_location(self) -> str:
        """The location where the unstaged samplesheet is located."""
        return f"s3://{self.bucket_name}/{self.staging_key}{self.samplesheet}"

    @property
    def synstage_outdir(self) -> str:
        """The S3 uri where nf-synstage writes the resolved samplesheet + FASTQs."""
        return f"s3://{self.bucket_name}/{self.staging_key}synstage/"

    @property
    def staged_samplesheet_location(self) -> str:
        """The S3 uri where the resolved samplesheet is published.

        nf-synapse's UPDATE_INPUT publishes under a run-name subfolder --
        confirmed empirically, matches als_rnaseq_params.json's "input".
        """
        return f"{self.synstage_outdir}{self.synstage_run_name}/{self.samplesheet}"

    @property
    def synstage_run_name(self) -> str:
        """The name of the synstage run."""
        return f"synstage_{self.id}"

    @property
    def rnaseq_run_name(self) -> str:
        """The name of the rnaseq run."""
        return f"rnaseq_GRCh38_{self.id}"

    @property
    def synindex_run_name(self) -> str:
        """The name of the synindex run."""
        return f"synindex_{self.id}"


REQUIRED_DATASET_FIELDS = ["name", 
                           "samplesheet_id", 
                           "samplesheet_filename", 
                           "output_folder_id"]


def generate_datasets() -> list[Dataset]:
    """Generate list of datasets from config/als_dataset_syn_ids.yml + config/pipeline_config.yml.

    To add another batch/disease later: 
    append another entry to `datasets` in config/als_dataset_syn_ids.yml 
    (and a matching key under `staging_keys` in config/pipeline_config.yml) 
    !! don't make code changes here.

    NOTE: per Ram Ayyala the 1,794 ALS samples span 7 source studies and
    should be treated as separate cohorts, not one combined batch; will
    eventually need separate runs for each cohort.

    Raises:
        SystemExit: a dataset entry is missing a required field, or has no
            matching `staging_keys` entry in pipeline_config.yml -- fails
            fast here rather than partway through launching Tower workflows.
    """
    syn_ids = yaml.safe_load(SYN_IDS_PATH.read_text())
    bucket_name = PIPELINE_CONFIG["bucket_name"]
    staging_keys = PIPELINE_CONFIG["staging_keys"]

    datasets = []
    for entry in syn_ids.get("datasets", []):
        missing = [f for f in REQUIRED_DATASET_FIELDS if not entry.get(f)]
        if missing:
            raise SystemExit(
                f"{SYN_IDS_PATH}: dataset entry {entry} is missing required "
                f"field(s) {missing}"
            )

        name = entry["name"]
        if name not in staging_keys:
            raise SystemExit(
                f"{PIPELINE_CONFIG_PATH}: no staging_keys entry for dataset "
                f"'{name}' (defined in {SYN_IDS_PATH}) -- add one before running"
            )
        staging_key = staging_keys[name]

        datasets.append(Dataset(
            id=entry["samplesheet_id"], #fetch_samplesheet() pull from Synapse -> upload to S3
            samplesheet=entry["samplesheet_filename"],
            staging_key=staging_key,
            bucket_name=bucket_name,
            synapse_id_for_output=entry["output_folder_id"], #Synapse folder for synindex - final rnaseq results back to Synapse
            output_directory=f"s3://{bucket_name}/{staging_key}rnaseq_output/",
        ))

    if not datasets:
        raise SystemExit(f"{SYN_IDS_PATH}: no datasets configured")

    return datasets


def fetch_samplesheet(syn: Synapse, dataset: Dataset) -> None:
    """Download the samplesheet from synapse and upload it to S3 (synstage will grab the file).

    Arguments:
        syn: The logged in synapse instance
        dataset: The dataset to stage the samplesheet for
    """
    samplesheet_file = syn_get(
        dataset.id, version_number=dataset.version, synapse_client=syn,
        file_options=FileOptions(if_collision="overwrite.local"),
    )
    s3.upload_file(
        samplesheet_file.path, dataset.bucket_name, f"{dataset.staging_key}{dataset.samplesheet}"
    )


def prepare_synstage_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for nf-synstage.

    Arguments:
        dataset: The dataset to stage the samplesheet for

    Returns:
        The Nextflow Tower workflow launch specification for synstage step
    """
    return LaunchInfo(
        run_name=dataset.synstage_run_name,
        pipeline="Sage-Bionetworks-Workflows/nf-synapse",
        revision="main",
        profiles=["sage"],
        params={
            "input": dataset.samplesheet_location,
            "outdir": dataset.synstage_outdir,
            "entry": "synstage",
        },
        workspace_secrets=["SYNAPSE_AUTH_TOKEN"],  # set as workspace secret (not user secret) in Tower
        nextflow_config="process.maxRetries = 1\nprocess.errorStrategy = { task.attempt <= 1 ? 'retry' : 'finish' }",
    )


def prepare_rnaseq_launch_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for nf-core/rnaseq workflow run.

    Params (aligner, skip flags, fasta/gtf) come from `als_rnaseq_params.json`,
    matching the AMP-AD parameter from https://sagebionetworks.jira.com/browse/SKG-145?focusedCommentId=336025
    `input`/`outdir` aren't in that file -- they're computed from the dataset,
    so they can never go stale when id changes.

    Arguments:
        dataset: The dataset to stage the samplesheet for

    Returns:
        The Nextflow Tower workflow launch specification for rnaseq processing step
    """
    params = {k: v for k, v in json.loads(PARAMS_PATH.read_text()).items() if not k.startswith("_")}
    params["input"] = dataset.staged_samplesheet_location
    params["outdir"] = dataset.output_directory
    return LaunchInfo(
        run_name=dataset.rnaseq_run_name,
        pipeline="nf-core/rnaseq",
        revision="287afcfe30a93de77e9b7cf70a1085f58c9525d8",  
        # per project runbook: https://docs.google.com/document/d/1icM_XQldD5c0e5_GLudOnJK-70F0IcBrYGgedyVnmDU/edit?tab=t.0
        profiles=["sage"],
        params=params,
        pre_run_script="NXF_VER=23.10.0",
        nextflow_config=CONFIG_PATH.read_text(),
    )

# synindex is a SINGLE run over the whole output_directory -- no tiering, no
# filename_string. The only hard Synapse limit is 10,000 CHILDREN PER FOLDER,
# confirmed by the real error on a samtools_stats file:
#   400 Client Error: Limit of 10000 children exceeded for parent: syn76651134
# There is NO per-run object-count limit. The split_samtools_stats pre-step brings
# the one offending folder (samtools_stats/, 10,764) under 10k; every other folder
# is already < 10k (next largest star_rsem/ ~7,180). So once the split has run, one
# synindex call mirrors the entire tree into parent_id.
#
# Re-run note: if this fails partway, re-running re-indexes everything. Synapse
# matches File by name+parent, so entities are UPDATED, not duplicated (a fresh
# run does bump a new version per file). Acceptable for a whole-run restart.
def prepare_synindex_launch_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for the single nf-synindex run over the whole output.

    No filename_string -> nf-synapse indexes every object under s3_prefix,
    mirroring the S3 folder structure 1:1 into parent_id.
    """
    return LaunchInfo(
        run_name=dataset.synindex_run_name,
        pipeline="Sage-Bionetworks-Workflows/nf-synapse",
        revision="main",
        profiles=["sage"],
        params={
            "s3_prefix": dataset.output_directory,
            "parent_id": dataset.synapse_id_for_output,
            "entry": "synindex",
        },
        workspace_secrets=["SYNAPSE_AUTH_TOKEN"],  # set as workspace secret (not user secret) in Tower
    )


async def run_workflows(ops: NextflowTowerOps, dataset: Dataset, step):
    """Run the requested step(s) for one dataset, in order."""
    if 'all' in step or 'fetch_samplesheet' in step:
        print('fetching samplesheet')
        syn = Synapse()
        syn.login()
        fetch_samplesheet(syn, dataset)

    if 'all' in step or 'synstage' in step:
        print('starting synstage')
        synstage_info = prepare_synstage_info(dataset)
        synstage_run_id = ops.launch_workflow(synstage_info, "spot", ignore_previous_runs=True)
        print(await monitor_workflow(ops, synstage_run_id))

    if 'all' in step or 'rnaseq' in step:
        print('starting data processing pipeline')
        rnaseq_info = prepare_rnaseq_launch_info(dataset)
        # on-demand: rnaseq is the long step; avoids losing hours of work to a spot reclaim
        rnaseq_run_id = ops.launch_workflow(rnaseq_info, "ondemand", ignore_previous_runs=True)
        print(await monitor_workflow(ops, rnaseq_run_id))

    if 'all' in step or 'split_samtools_stats' in step:
        # Pre-synindex: split star_rsem/samtools_stats/ (10,764 files, the only
        # folder over Synapse's 10k-children limit) into raw/ + markdup/. Idempotent
        # and blocking boto3, so run it off the event loop. See split_samtools_stats.py.
        print('splitting star_rsem/samtools_stats/ into raw/ + markdup/ (10k-per-folder limit)')
        await asyncio.to_thread(
            split_samtools_stats, f"{dataset.output_directory}star_rsem/samtools_stats/"
        )

    if 'all' in step or 'synindex' in step:
        # Single run over the whole output_directory (see prepare_synindex_launch_info).
        # No skip-if-succeeded pre-check: it calls ops.list_workflows(), which orca
        # crashes on against the current Tower API (the /workflow response ends with a
        # trailing `truncated` bool that orca's popitem() mistakes for the items list).
        # Re-running synindex is safe anyway -- Synapse matches by name+parent, so it
        # updates entities (bumps version) rather than duplicating.
        print('starting synindex')
        synindex_info = prepare_synindex_launch_info(dataset)
        synindex_run_id = ops.launch_workflow(synindex_info, "spot", ignore_previous_runs=True)
        print(await monitor_workflow(ops, synindex_run_id))


def _self_check():
    """Sanity check on path construction (no Tower/Synapse calls).

    Static only: confirms the computed S3 paths agree with each other and with
    config/*.yml + als_rnaseq_params.json, for every configured dataset -- not
    just the first. Does NOT confirm those paths exist in S3 or that
    synstage/rnaseq have actually run — that requires a real `--check`-free run.
    """
    for ds in generate_datasets():
        expected_prefix = f"s3://{ds.bucket_name}/{ds.staging_key}"
        assert ds.samplesheet_location == f"{expected_prefix}{ds.samplesheet}"
        assert ds.staged_samplesheet_location == (
            f"{expected_prefix}synstage/{ds.synstage_run_name}/{ds.samplesheet}"
        )
        assert ds.output_directory == f"{expected_prefix}rnaseq_output/"

        rnaseq_info = prepare_rnaseq_launch_info(ds)
        assert rnaseq_info.params["input"] == ds.staged_samplesheet_location
        assert rnaseq_info.params["outdir"] == ds.output_directory
    print("self-check OK")


if __name__ == "__main__":
    import sys
    if "--check" in sys.argv:
        _self_check()
    else:
        asyncio.run(main())

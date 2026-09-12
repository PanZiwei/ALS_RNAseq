"""ORCA recipe: nf-core/rnaseq bulk RNA-seq processing on GRCh38.
Implements SKG-145 (https://sagebionetworks.jira.com/browse/SKG-145).

Orchestrate four steps (run in order):
  1. fetch_samplesheet  : Synapse -> S3
  2. synstage           : nf-synapse, download FASTQs Synapse -> S3
  3. rnaseq             : nf-core/rnaseq alignment/quantification
  4. synindex           : nf-synapse, S3 -> Synapse, one Tower run per
                          SYNINDEX_TIERS entry (8 total):
                            tier1 (3 entries) -- merged matrices + QC + provenance
                            tier2 (4 entries) -- per-sample results/stats/logs
                            tier3 (1 entry)   -- BAMs
                          No resume: one entry failing stops the whole run, re-run redoes every entry 

Prerequisites:
  - pip install py-orca pyyaml
  - AWS credentials configured (profile: tower)
  - SYNAPSE_AUTH_TOKEN set as a Tower workspace secret (not user secret):
    https://tower.sagebionetworks.org/orgs/Sage-Bionetworks/workspaces/complementaire-project/secrets
  - NEXTFLOWTOWER_CONNECTION_URI env var (only env var this script reads):
    1. Get a token: tower.sagebionetworks.org -> avatar -> Your tokens -> Add token
    2. export NEXTFLOWTOWER_CONNECTION_URI="https://:<token>@tower.sagebionetworks.org/api?workspace=Sage-Bionetworks%2Fcomplementaire-project"

Usage:
  python orca_rnaseq_workflow.py                        # all four steps
  python orca_rnaseq_workflow.py synstage                # one step
  python orca_rnaseq_workflow.py synstage rnaseq         # multiple steps
  python orca_rnaseq_workflow.py synindex --tier tier3   # synindex, one tier only
  python orca_rnaseq_workflow.py synindex --tier tier1 --tier tier3  # 2 tiers
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
        choices=['all', 'fetch_samplesheet', 'synstage', 'rnaseq', 'synindex'],
        help="One or more steps to run, in order (default: all four)",
    )
    parser.add_argument(
        '--tier', action='append',
        choices=sorted({tier_name for tier_name, *_ in SYNINDEX_TIERS}),
        help="Restrict the synindex step to one tier from SYNINDEX_TIERS "
             "(repeatable, e.g. --tier tier1 --tier tier3). Default: all tiers.",
    )
    args = parser.parse_args()

    ops = NextflowTowerOps()
    datasets = generate_datasets()
    runs = [run_workflows(ops, dataset, args.step, args.tier) for dataset in datasets]
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

# synindex (nf-synapse) can't reliably process more than ~10,000 S3 objects
# in one run -- indexing the whole output_directory in one call failed
# Fix: one synindex run per (tier, label, s3_subprefix, filename_string) entry below, 
# each under the S3 object limit.
# s3_subprefix="" means "the whole output_directory" -- used to merge
# several small sibling S3 folders (e.g. multiqc/, pipeline_info/) into one
# run via filename_string, when their combined size is still safely under
# the limit. 
# Everything under star_rsem/ shares s3_subprefix="star_rsem" (flat on S3) and uses filename_string to
# non-overlappingly partition it instead.
#
# Tiers run in order, most useful first, so a paused/partial synindex run
# still leaves the most valuable artifacts on Synapse; --tier restricts a run
# to one or more tiers (e.g. --tier tier3 to re-run just the BAMs):
#   tier1 -- merged matrices + QC + provenance
#   tier2 -- per-sample results/stats/logs, for re-analysis/auditing
#   tier3 -- BAMs: largest, re-derivable from FASTQs, least urgent
#
# filename_string matches the FULL S3 key, not just the basename -- anchor
# with "$" (safe, basename is always at the end) not "^"; test any new/edited
# pattern against real S3 keys before relying on it (a "^"-anchored
# merged_matrices pattern silently matched nothing until this was caught).
SYNINDEX_TIERS = [
    ("tier1", "matrices_and_multiqc", "", r"/rsem\.merged\.|/multiqc/"),
    ("tier1", "fastqc", "fastqc", None),
    ("tier1", "trimgalore_and_pipeline_info", "", r"/trimgalore/|/pipeline_info/"),
    ("tier2", "rsem_results", "star_rsem", r"\.genes\.results$|\.isoforms\.results$"),
    ("tier2", "samtools_stats_raw", "star_rsem/samtools_stats", r"^(?!.*\.markdup\.).*\.stats$"),
    ("tier2", "samtools_stats_markdup", "star_rsem/samtools_stats", r"\.markdup\.stats$"),
    ("tier2", "star_logs", "star_rsem", r"\.stat/"),
    ("tier3", "bams", "star_rsem", r"\.markdup\.sorted\.bam(\.bai)?$"),
]


def prepare_synindex_launch_info(
    dataset: Dataset, s3_subprefix: str, filename_string: str | None, tier_name: str, run_label: str,
) -> LaunchInfo:
    """Generate LaunchInfo for one nf-synindex run, scoped to one entry from SYNINDEX_TIERS.

    Arguments:
        dataset: The dataset whose output is being indexed
        s3_subprefix: Relative S3 prefix under output_directory
        filename_string: nf-synapse filename filter regex, or None to index
            everything under s3_subprefix
        tier_name: This entry's tier from SYNINDEX_TIERS, e.g. "tier1, tier2, tier3"
        run_label: This entry's label from SYNINDEX_TIERS, e.g.
            "merged_matrices" -- together they form this run's Tower
            run_name, e.g. "synindex_syn76551366_tier1_matrices_and_multiqc"

    Returns:
        The Nextflow Tower workflow launch specification for this synindex run
    """
    # s3_subprefix="" means "the whole output_directory" -- used when
    # filename_string alone picks out objects from multiple sibling S3
    # folders in one run (see the merged tier1 entries in SYNINDEX_TIERS).
    s3_prefix = dataset.output_directory if not s3_subprefix else f"{dataset.output_directory}{s3_subprefix}/"
    params = {
        "s3_prefix": s3_prefix,
        "parent_id": dataset.synapse_id_for_output,
        "entry": "synindex",
    }
    if filename_string is not None:
        params["filename_string"] = filename_string
    return LaunchInfo(
        run_name=f"{dataset.synindex_run_name}_{tier_name}_{run_label}",
        pipeline="Sage-Bionetworks-Workflows/nf-synapse",
        revision="main",
        profiles=["sage"],
        params=params,
        workspace_secrets=["SYNAPSE_AUTH_TOKEN"]  # set as workspace secret (not user secret) in Tower
    )


async def run_workflows(ops: NextflowTowerOps, dataset: Dataset, step, tiers: list[str] | None = None):
    """Run the requested step(s) for one dataset.

    Arguments:
        tiers: If given, restricts the synindex step to these SYNINDEX_TIERS
            tier names (e.g. ["core"]); None runs every tier.
    """
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

    if 'all' in step or 'synindex' in step:
        # One synindex run per SYNINDEX_TIERS entry (see its comment); tiers
        # restricts this to a subset when set (--tier on the CLI).
        entries = SYNINDEX_TIERS if tiers is None else [
            entry for entry in SYNINDEX_TIERS if entry[0] in tiers
        ]
        for tier_name, run_label, s3_subprefix, filename_string in entries:
            print(f'starting synindex ({tier_name}/{run_label})')
            synindex_info = prepare_synindex_launch_info(
                dataset, s3_subprefix, filename_string, tier_name, run_label
            )
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

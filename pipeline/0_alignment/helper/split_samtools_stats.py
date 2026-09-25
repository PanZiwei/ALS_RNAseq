"""Split star_rsem/samtools_stats/ into raw/ and markdup/ subfolders.

Pre-synindex step. samtools_stats/ holds 6 files/sample, so at 1,794 samples it
has 10,764 objects in one folder -- over Synapse's hard 10,000-children-per-folder
limit, which fails synindex on that folder. It is the ONLY output folder over the
limit (next largest, star_rsem/, is ~7,180). Splitting by type halves it to 5,382
each, both safely under the limit, and keeps every file individually browsable on
Synapse (vs. zipping into one archive).

Uses boto3 per-object copy/delete (server-side, no download). Idempotent: a file
already at its correct destination is skipped, and the whole step no-ops if the
folder is already under the limit.

Config (bucket, AWS profile) is read from config/aws_s3_config.yml, same as
orca_rnaseq_workflow.py -- not hardcoded.

Usage:
  python split_samtools_stats.py --check                 # static self-check, no AWS
  python split_samtools_stats.py --dry-run               # list what would move, no writes
  python split_samtools_stats.py                         # split the default run's output
  python split_samtools_stats.py --s3-prefix s3://bucket/path/star_rsem/samtools_stats/
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
import yaml

CONFIG_DIR = Path(__file__).parent.parent.parent.parent / "config"
AWS_S3_CONFIG = yaml.safe_load((CONFIG_DIR / "aws_s3_config.yml").read_text())

AWS_PROFILE = AWS_S3_CONFIG["aws_profile"]
BUCKET = AWS_S3_CONFIG["bucket_name"]
# Default target: the ALS combined run's samtools_stats/. --s3-prefix overrides.
DEFAULT_S3_PREFIX = f"s3://{BUCKET}/ALS/rnaseq_output/star_rsem/samtools_stats/"

MAX_WORKERS = 20
PROGRESS_EVERY = 500

session = boto3.Session(profile_name=AWS_PROFILE)
s3 = session.client("s3")


def _parse_s3_prefix(s3_prefix: str) -> tuple[str, str]:
    bucket, _, key_prefix = s3_prefix.removeprefix("s3://").partition("/")
    return bucket, key_prefix


def _list_all_objects(bucket: str, prefix: str) -> list[str]:
    keys = []
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            keys.append(obj["Key"])
    return keys


def _target_key(key: str, samtools_stats_prefix: str) -> str:
    """raw/ vs markdup/ by basename. markdup files carry '.markdup.sorted.bam.';
    using the basename (not the full key) makes this idempotent regardless of any
    prior subfolder nesting."""
    basename = key.rsplit("/", 1)[-1]
    subfolder = "markdup" if ".markdup.sorted.bam." in basename else "raw"
    return f"{samtools_stats_prefix}{subfolder}/{basename}"


def _move_object(bucket: str, key: str, target_key: str) -> None:
    s3.copy_object(Bucket=bucket, CopySource={"Bucket": bucket, "Key": key}, Key=target_key)
    s3.delete_object(Bucket=bucket, Key=key)


def split_samtools_stats(samtools_stats_prefix: str, dry_run: bool = False) -> None:
    """Arguments: samtools_stats_prefix -- S3 URI of the samtools_stats/ folder itself."""
    bucket, key_prefix = _parse_s3_prefix(samtools_stats_prefix)

    all_keys = _list_all_objects(bucket, key_prefix)
    files = [k for k in all_keys if not k.endswith("/")]
    to_move = [(k, _target_key(k, key_prefix)) for k in files if k != _target_key(k, key_prefix)]
    print(f"samtools_stats/: {len(files)} files, {len(to_move)} need moving (rest already split)")

    if dry_run:
        for k, target in to_move[:5]:
            print(f"  would move {k} -> {target}")
        print("dry-run only -- no changes made")
        return
    if not to_move:
        print("Nothing to move -- already split.")
        return

    done = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(_move_object, bucket, k, target): k for k, target in to_move}
        for future in as_completed(futures):
            future.result()  # raises if the move failed -- fail loud
            done += 1
            if done % PROGRESS_EVERY == 0 or done == len(to_move):
                print(f"  {done}/{len(to_move)} moved")

    after = _list_all_objects(bucket, key_prefix)
    raw = sum(1 for k in after if f"{key_prefix}raw/" in k and not k.endswith("/"))
    markdup = sum(1 for k in after if f"{key_prefix}markdup/" in k and not k.endswith("/"))
    print(f"raw/: {raw} files, markdup/: {markdup} files")
    assert raw < 10000 and markdup < 10000, "still over the 10k-per-folder limit -- investigate"
    print("split OK")


def _self_check():
    """Static check only -- confirms path parsing + raw/markdup routing, no S3 calls."""
    bucket, key_prefix = _parse_s3_prefix(DEFAULT_S3_PREFIX)
    assert bucket == BUCKET
    assert _target_key(f"{key_prefix}SRR1.sorted.bam.stats", key_prefix) == \
        f"{key_prefix}raw/SRR1.sorted.bam.stats"
    assert _target_key(f"{key_prefix}SRR1.markdup.sorted.bam.stats", key_prefix) == \
        f"{key_prefix}markdup/SRR1.markdup.sorted.bam.stats"
    # idempotent regardless of prior nesting depth -- basename is all that matters
    assert _target_key(f"{key_prefix}markdup/SRR1.markdup.sorted.bam.stats", key_prefix) == \
        f"{key_prefix}markdup/SRR1.markdup.sorted.bam.stats"
    print("self-check OK")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--s3-prefix", default=DEFAULT_S3_PREFIX, help="S3 URI of samtools_stats/ itself")
    parser.add_argument("--check", action="store_true", help="Static self-check only, no AWS calls")
    parser.add_argument("--dry-run", action="store_true", help="List what would move, make no changes")
    args = parser.parse_args()

    if args.check:
        _self_check()
    else:
        split_samtools_stats(args.s3_prefix, dry_run=args.dry_run)

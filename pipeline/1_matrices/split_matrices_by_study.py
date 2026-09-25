"""Split the combined run's merged RSEM matrices into one file per study.

The split streams each TSV line by line and writes all study files in one pass,
so the ~2.2 GiB matrices are never loaded into RAM.

Runs on EC2 (controlled-access data). Set a Synapse token first -- a personal
access token with Download scope (add Modify for `upload`):
  export SYNAPSE_AUTH_TOKEN='<token>'

Steps (run one, several in order, or default `split`):
  fetch   download all inputs into --data-dir: the 4 merged matrices + the two
          sample->study mapping files (dataset_targets, samplesheet_with_dataset_names)
  split   split each matrix in --data-dir into one file per study, into --out-dir
  upload  upload the per-study files to Synapse (--upload-parent synXXXXXXX)

Usage:
  python split_matrices_by_study.py --check    # verify mapping only
  python split_matrices_by_study.py fetch      # download the 4 matrices
  python split_matrices_by_study.py split
  python split_matrices_by_study.py upload --upload-parent synXXXXXXX

All source Synapse IDs (matrices folder + the two mapping CSVs) live in
config/als_dataset_syn_ids.yml's split_matrices_inputs. The mapping CSVs are
pulled automatically by split/upload/--check; the matrices by `fetch`. All are
cached in --data-dir.
"""
import argparse
import csv
import sys
from operator import itemgetter
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).parent.parent.parent
CONFIG = yaml.safe_load((REPO_ROOT / "config" / "als_dataset_syn_ids.yml").read_text())
SPLIT_INPUTS = CONFIG["split_matrices_inputs"]

DEFAULT_DATA_DIR = REPO_ROOT / "data" / "rsem_merged_tsv"
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "per_study"

MATRIX_FILES = [
    "rsem.merged.gene_counts.tsv",
    "rsem.merged.gene_tpm.tsv",
    "rsem.merged.transcript_counts.tsv",
    "rsem.merged.transcript_tpm.tsv",
]
EXPECTED_SAMPLES = 1794
STUDY_KEYWORDS = SPLIT_INPUTS["study_keywords"]  # dataset_synid -> filename slug


def _synapse():
    import synapseclient
    syn = synapseclient.Synapse()
    syn.login(silent=True)
    return syn


def _children(syn, parent_id: str) -> list[dict]:
    """List a folder's children via the non-deprecated (async) get_children API,
    run synchronously so the rest of the script stays plain."""
    import asyncio
    from synapseclient.api import get_children

    async def _collect():
        return [c async for c in get_children(parent_id, synapse_client=syn)]

    return asyncio.run(_collect())


def _fetch(synapse_id: str, filename: str, data_dir: Path) -> Path:
    """Download a Synapse file into data_dir as `filename`, cached there.
    Streams to disk (downloadLocation) so multi-GB matrices never load into RAM."""
    from synapseclient.operations import FileOptions
    from synapseclient.operations import get as syn_get

    local = data_dir / filename
    if local.exists():
        return local
    data_dir.mkdir(parents=True, exist_ok=True)
    entity = syn_get(
        synapse_id,
        file_options=FileOptions(download_location=str(data_dir)),
        synapse_client=_synapse(),
    )
    got = Path(entity.path)
    if got != local:
        got.replace(local)
    return local


def _load_studies(data_dir: Path) -> dict[str, str]:
    """dataset_synid -> dataset_name, validated against STUDY_KEYWORDS."""
    path = _fetch(SPLIT_INPUTS["dataset_targets_id"], "dataset_targets.csv", data_dir)
    with open(path) as f:
        studies = {r["dataset_synid"]: r["dataset_name"] for r in csv.DictReader(f)}
    if set(studies) != set(STUDY_KEYWORDS):
        raise SystemExit(
            f"dataset_targets.csv study IDs {sorted(studies)} do not match "
            f"STUDY_KEYWORDS {sorted(STUDY_KEYWORDS)} -- update config"
        )
    return studies


def build_sample_to_study(data_dir: Path) -> dict[str, str]:
    """Map each sample -> study_id, joining the samplesheet on dataset_name."""
    name_to_synid = {name: sid for sid, name in _load_studies(data_dir).items()}
    samplesheet = _fetch(
        SPLIT_INPUTS["samplesheet_with_dataset_names_id"],
        "samplesheet_with_dataset_names.csv", data_dir,
    )
    with open(samplesheet) as f:
        return {r["sample"]: name_to_synid[r["dataset"]] for r in csv.DictReader(f)}


def out_name(study_id: str, matrix_name: str) -> str:
    return f"{study_id}_{STUDY_KEYWORDS[study_id]}.{matrix_name}"


def split_matrix(src: Path, sample_to_study: dict[str, str], out_dir: Path) -> None:
    """Stream one merged matrix into one file per study (single pass, low memory)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(src) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        unknown = [c for c in header[2:] if c not in sample_to_study]
        if unknown:
            raise SystemExit(f"{src.name}: {len(unknown)} column(s) not in samplesheet, e.g. {unknown[:5]}")

        # Per study: the id columns (0,1) + that study's sample columns, picked
        # with itemgetter (C-speed) so the row loop has no per-column Python work.
        writers = []  # (handle, itemgetter)
        for study_id in STUDY_KEYWORDS:
            cols = [0, 1] + [i for i, c in enumerate(header) if i >= 2 and sample_to_study[c] == study_id]
            if len(cols) == 2:
                continue
            fh_out = open(out_dir / out_name(study_id, src.name), "w")
            pick = itemgetter(*cols)
            fh_out.write("\t".join(pick(header)) + "\n")
            writers.append((fh_out, pick))
        try:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                for fh_out, pick in writers:
                    fh_out.write("\t".join(pick(parts)) + "\n")
        finally:
            for fh_out, _ in writers:
                fh_out.close()
    print(f"  {src.name}: {len(writers)} study files")


# -- steps --------------------------------------------------------------------

def cmd_fetch(args) -> None:
    """Download all inputs into --data-dir: the 4 merged matrices (matched by name
    from the source folder, which holds other files too) plus the two sample->study
    mapping files (dataset_targets.csv, samplesheet_with_dataset_names.csv). Folder
    ID defaults to split_matrices_inputs.rsem_matrices_folder_id; --synid overrides."""
    synid = args.synid or SPLIT_INPUTS.get("rsem_matrices_folder_id")
    if not synid:
        raise SystemExit("no folder ID -- set rsem_matrices_folder_id in config or pass --synid")
    ids = {c["name"]: c["id"] for c in _children(_synapse(), synid)}
    for name in MATRIX_FILES:
        if name not in ids:
            raise SystemExit(f"{name} not found in {synid}")
        print(f"fetching {name} ({ids[name]}) ...")
        _fetch(ids[name], name, args.data_dir)
    # also grab the sample->study mapping files, so all inputs are local after `fetch`
    _fetch(SPLIT_INPUTS["dataset_targets_id"], "dataset_targets.csv", args.data_dir)
    _fetch(SPLIT_INPUTS["samplesheet_with_dataset_names_id"],
           "samplesheet_with_dataset_names.csv", args.data_dir)
    print("fetched 4 merged matrices + mapping files "
          "(dataset_targets.csv, samplesheet_with_dataset_names.csv)")


def cmd_split(args) -> None:
    sample_to_study = build_sample_to_study(args.data_dir)
    for filename in MATRIX_FILES:
        src = args.data_dir / filename
        if not src.exists():
            raise SystemExit(f"{src} not found -- run `fetch --synid ...` first")
        split_matrix(src, sample_to_study, args.out_dir)


def cmd_upload(args) -> None:
    """Upload per-study matrices to Synapse: one folder per study, each file
    annotated with its dataset. Re-runs are safe -- store() matches by
    name+parent, so it updates entities (new version) rather than duplicating."""
    if not args.upload_parent:
        raise SystemExit("upload requires --upload-parent synXXXXXXX")
    from synapseclient import File, Folder
    syn = _synapse()
    for study_id, dataset_name in _load_studies(args.data_dir).items():
        paths = [args.out_dir / out_name(study_id, m) for m in MATRIX_FILES]
        missing = [p for p in paths if not p.exists()]
        if missing:
            raise SystemExit(f"{[str(p) for p in missing]} not found -- run `split` first")
        folder_name = f"{study_id}_{STUDY_KEYWORDS[study_id]}"
        print(f"{folder_name} ({dataset_name[:60]}...)")
        if args.dry_run:
            for p in paths:
                print(f"  would upload {p.name} -> {args.upload_parent}/{folder_name}")
            continue
        folder = syn.store(Folder(folder_name, parent=args.upload_parent))
        for p in paths:
            e = syn.store(File(str(p), parent=folder,
                               annotations={"dataset_synid": study_id, "dataset_name": dataset_name}))
            print(f"  uploaded {p.name} -> {e.id}")


def _self_check(data_dir: Path) -> None:
    """Verify the sample->study mapping resolves to the expected study set/count."""
    studies = _load_studies(data_dir)
    sample_to_study = build_sample_to_study(data_dir)
    counts: dict[str, int] = {}
    for sid in sample_to_study.values():
        counts[sid] = counts.get(sid, 0) + 1
    assert len(sample_to_study) == EXPECTED_SAMPLES, f"expected {EXPECTED_SAMPLES}, got {len(sample_to_study)}"
    assert set(counts) == set(STUDY_KEYWORDS), f"study mismatch: {set(counts)} vs {set(STUDY_KEYWORDS)}"
    for sid, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"{sid} ({STUDY_KEYWORDS[sid]}): {n} samples -- {studies[sid][:60]}")
    print("self-check OK")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("step", nargs="*", help="fetch, split, upload (default: split)")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help=f"merged matrices + cached CSVs (default: {DEFAULT_DATA_DIR})")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help=f"per-study output (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--synid", help="fetch: override the config's rsem_matrices_folder_id")
    parser.add_argument("--upload-parent", help="upload: Synapse folder to create per-study folders under")
    parser.add_argument("--dry-run", action="store_true", help="upload: print, store nothing")
    parser.add_argument("--check", action="store_true", help="verify the sample->study mapping and exit")
    args = parser.parse_args()

    if args.check:
        _self_check(args.data_dir)
        return

    steps = {"fetch": cmd_fetch, "split": cmd_split, "upload": cmd_upload}
    for name in args.step or ["split"]:
        if name not in steps:
            raise SystemExit(f"unknown step {name!r} (choose from: fetch, split, upload)")
        print(f"== {name} ==")
        steps[name](args)


if __name__ == "__main__":
    sys.exit(main())

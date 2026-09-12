"""Split the combined run's merged RSEM matrices into one file per study.

The split streams the TSV line by line and writes all 7 study files in a
single pass. rsem.merged.transcript_*.tsv is ~2.2 GiB; loading that into
pandas would need well over 10 GiB of RAM, more than this host has.

Usage:
  python split_merged_matrics_by_study.py --check      # verify sample->study mapping, no I/O
  python split_merged_matrics_by_study.py --test       # split 1 matrix, verify 2 real samples
                                                       # match the source values exactly
  python split_merged_matrics_by_study.py split        # split all 4 matrices x 7 studies
  python split_merged_matrics_by_study.py upload --upload-parent synXXXXXXX
                                                       # upload per-study files to Synapse
  python split_merged_matrics_by_study.py split upload --upload-parent synXXXXXXX

The merged matrices are expected to be already downloaded into --data-dir.

samplesheet_with_dataset_names.csv and dataset_targets.csv are not committed
to this repo -- Synapse is their source of truth (see
config/als_dataset_syn_ids.yml's split_matrices_inputs). They're downloaded
into --data-dir on first use and cached there; delete the cached copy to
force a re-download if the Synapse files change.
"""
import argparse
import csv
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).parent.parent.parent
CONFIG = yaml.safe_load((REPO_ROOT / "config" / "als_dataset_syn_ids.yml").read_text())
SPLIT_INPUTS = CONFIG["split_matrices_inputs"]

DEFAULT_DATA_DIR = REPO_ROOT / "data" / "rsem_merged_tsv"
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "per_study"


def _fetch_from_synapse(synapse_id: str, filename: str, data_dir: Path) -> Path:
    """Download a file from Synapse into data_dir, caching it there.

    Arguments:
        synapse_id: The Synapse ID of the file to download
        filename: Local filename to cache it as, under data_dir
        data_dir: Directory to cache the downloaded file in

    Returns:
        Path to the (possibly cached) local file
    """
    local_path = data_dir / filename
    if local_path.exists():
        return local_path

    import synapseclient
    from synapseclient.operations import get as syn_get

    data_dir.mkdir(parents=True, exist_ok=True)
    syn = synapseclient.Synapse()
    syn.login(silent=True)
    downloaded = syn_get(synapse_id, synapse_client=syn)
    local_path.write_bytes(Path(downloaded.path).read_bytes())
    return local_path

MATRIX_FILES = [
    "rsem.merged.gene_counts.tsv",
    "rsem.merged.gene_tpm.tsv",
    "rsem.merged.transcript_counts.tsv",
    "rsem.merged.transcript_tpm.tsv",
]

EXPECTED_SAMPLES = 1794

# Short filename slug per study -- see config/als_dataset_syn_ids.yml's
# split_matrices_inputs.study_keywords for the values and rationale.
STUDY_KEYWORDS = SPLIT_INPUTS["study_keywords"]


def _load_studies(data_dir: Path) -> dict[str, str]:
    """dataset_synid -> dataset_name, validated against STUDY_KEYWORDS."""
    dataset_targets = _fetch_from_synapse(
        SPLIT_INPUTS["dataset_targets_id"], "dataset_targets.csv", data_dir
    )
    studies = {}
    with open(dataset_targets) as f:
        for row in csv.DictReader(f):
            studies[row["dataset_synid"]] = row["dataset_name"]
    if set(studies) != set(STUDY_KEYWORDS):
        raise SystemExit(
            f"dataset_targets.csv study IDs {sorted(studies)} do not match "
            f"STUDY_KEYWORDS {sorted(STUDY_KEYWORDS)} -- update STUDY_KEYWORDS"
        )
    return studies


def build_sample_to_study(data_dir: Path) -> dict[str, str]:
    """Map each sample to its study_id, via dataset_name as the join key."""
    dataset_targets = _fetch_from_synapse(
        SPLIT_INPUTS["dataset_targets_id"], "dataset_targets.csv", data_dir
    )
    samplesheet_with_dataset = _fetch_from_synapse(
        SPLIT_INPUTS["samplesheet_with_dataset_names_id"],
        "samplesheet_with_dataset_names.csv", data_dir,
    )

    name_to_synid = {}
    with open(dataset_targets) as f:
        for row in csv.DictReader(f):
            name_to_synid[row["dataset_name"]] = row["dataset_synid"]

    sample_to_study = {}
    with open(samplesheet_with_dataset) as f:
        for row in csv.DictReader(f):
            study_id = name_to_synid[row["dataset"]]
            sample_to_study[row["sample"]] = study_id
    return sample_to_study


def out_name(study_id: str, matrix_name: str) -> str:
    return f"{study_id}_{STUDY_KEYWORDS[study_id]}.{matrix_name}"


def split_matrix(local_path: Path, sample_to_study: dict[str, str], out_dir: Path) -> list[Path]:
    """Stream one merged matrix into one file per study. Returns written paths.

    Single pass: every study's output file is open at once and each input line
    is fanned out to all of them, so the 2.2 GiB matrices are never held in
    memory and are read from disk only once.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(local_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        id_cols = header[:2]  # gene_id/transcript_id + second RSEM column
        sample_cols = header[2:]

        unknown = [c for c in sample_cols if c not in sample_to_study]
        if unknown:
            raise SystemExit(
                f"{local_path.name}: {len(unknown)} column(s) not in the samplesheet, "
                f"e.g. {unknown[:5]}"
            )

        # Column indices per study, so the row loop is a plain list lookup.
        indices: dict[str, list[int]] = {}
        for study_id in STUDY_KEYWORDS:
            cols = [i for i, c in enumerate(header) if i >= 2 and sample_to_study[c] == study_id]
            if cols:
                indices[study_id] = [0, 1] + cols

        handles = {}
        written = []
        try:
            for study_id, idx in indices.items():
                path = out_dir / out_name(study_id, local_path.name)
                handles[study_id] = open(path, "w")
                handles[study_id].write("\t".join(header[i] for i in idx) + "\n")
                written.append(path)

            for line in fh:
                parts = line.rstrip("\n").split("\t")
                for study_id, idx in indices.items():
                    handles[study_id].write("\t".join(parts[i] for i in idx) + "\n")
        finally:
            for h in handles.values():
                h.close()

    for study_id, idx in indices.items():
        print(f"  {out_name(study_id, local_path.name)}: {len(idx) - len(id_cols)} samples")
    return written


# -- steps --------------------------------------------------------------------

def cmd_split(args) -> None:
    sample_to_study = build_sample_to_study(args.data_dir)
    for filename in MATRIX_FILES:
        local_path = args.data_dir / filename
        if not local_path.exists():
            raise SystemExit(f"{local_path} not found -- put the merged matrices in --data-dir")
        print(f"{filename}:")
        split_matrix(local_path, sample_to_study, args.out_dir)


def cmd_upload(args) -> None:
    """Upload the per-study matrices to Synapse, one folder per study."""
    if not args.upload_parent:
        raise SystemExit("upload requires --upload-parent synXXXXXXX")

    import synapseclient
    from synapseclient import File, Folder

    studies = _load_studies(args.data_dir)
    syn = synapseclient.Synapse()
    syn.login(silent=True)

    for study_id, dataset_name in studies.items():
        paths = [args.out_dir / out_name(study_id, m) for m in MATRIX_FILES]
        missing = [p for p in paths if not p.exists()]
        if missing:
            raise SystemExit(f"{[str(p) for p in missing]} not found -- run the `split` step first")

        folder_name = f"{study_id}_{STUDY_KEYWORDS[study_id]}"
        print(f"{folder_name} ({dataset_name[:60]}...)")
        if args.dry_run:
            for path in paths:
                print(f"  would upload {path.name} -> {args.upload_parent}/{folder_name}")
            continue

        folder = syn.store(Folder(folder_name, parent=args.upload_parent))
        for path in paths:
            entity = syn.store(
                File(
                    str(path),
                    parent=folder,
                    annotations={"dataset_synid": study_id, "dataset_name": dataset_name},
                )
            )
            print(f"  uploaded {path.name} -> {entity.id}")


# -- checks -------------------------------------------------------------------

def _self_check(data_dir: Path) -> None:
    studies = _load_studies(data_dir)
    sample_to_study = build_sample_to_study(data_dir)
    assert len(sample_to_study) == EXPECTED_SAMPLES, (
        f"expected {EXPECTED_SAMPLES} samples, got {len(sample_to_study)}"
    )
    counts: dict[str, int] = {}
    for study_id in sample_to_study.values():
        counts[study_id] = counts.get(study_id, 0) + 1
    assert sum(counts.values()) == EXPECTED_SAMPLES
    assert set(counts) == set(STUDY_KEYWORDS), (
        f"study_id mismatch: {set(counts)} vs {set(STUDY_KEYWORDS)}"
    )
    for study_id, count in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"{study_id} ({STUDY_KEYWORDS[study_id]}): {count} samples -- {studies[study_id][:60]}")
    print("self-check OK")


def _end_to_end_test(args) -> None:
    """Real-data correctness test: split 1 matrix, then verify the per-study
    output is identical to the source for 2 real sample columns -- not just
    structurally shaped right, but actually correct values."""
    from collections import Counter
    import pandas as pd

    sample_to_study = build_sample_to_study(args.data_dir)
    out_dir = Path("/tmp/merged_matrices_test")

    filename = "rsem.merged.gene_counts.tsv"  # smallest of the 4
    local_path = args.data_dir / filename
    if not local_path.exists():
        raise SystemExit(f"{local_path} not found -- put the merged matrices in --data-dir")

    df = pd.read_csv(local_path, sep="\t")
    id_cols = list(df.columns[:2])

    # pick 2 real samples from the same study, so the test also confirms
    # per-study grouping works, not just column selection
    study_id, _ = Counter(sample_to_study.values()).most_common(1)[0]
    test_samples = [s for s, st in sample_to_study.items() if st == study_id
                    and s in df.columns][:2]
    assert len(test_samples) == 2, f"couldn't find 2 test samples for {study_id}"
    print(f"testing samples: {test_samples} (study {study_id})")

    split_matrix(local_path, sample_to_study, out_dir)
    split_df = pd.read_csv(out_dir / out_name(study_id, filename), sep="\t")

    for sample in test_samples:
        original_col = df.set_index(id_cols[0])[sample]
        split_col = split_df.set_index(id_cols[0])[sample]
        pd.testing.assert_series_equal(original_col, split_col, check_names=False)
        print(f"  {sample}: values match exactly ({len(original_col)} genes)")

    print("end-to-end test OK -- split output matches source for tested samples")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # default=None, not a list: argparse validates a list default against
    # `choices` as a single value and rejects it.
    parser.add_argument("step", nargs="*", choices=["split", "upload"],
                        help="Steps to run, in order (default: split)")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                        help=f"Where the merged matrices live (default: {DEFAULT_DATA_DIR})")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help=f"Where per-study matrices are written (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--upload-parent",
                        help="Synapse ID of the folder to create per-study folders under")
    parser.add_argument("--dry-run", action="store_true",
                        help="For upload: print what would be uploaded, store nothing")
    parser.add_argument("--check", action="store_true",
                        help="Verify the sample->study mapping and exit")
    parser.add_argument("--test", action="store_true",
                        help="Run the end-to-end correctness test and exit")
    args = parser.parse_args()

    if args.check:
        _self_check(args.data_dir)
        return
    if args.test:
        _end_to_end_test(args)
        return

    steps = {"split": cmd_split, "upload": cmd_upload}
    for name in args.step or ["split"]:
        print(f"== {name} ==")
        steps[name](args)


if __name__ == "__main__":
    sys.exit(main())

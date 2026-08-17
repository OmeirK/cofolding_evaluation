#!/usr/bin/env python3
"""Prepare pose-level co-folding results for downstream plotting.

This script combines one or more previously analysed co-folding result tables,
adds annotations from the EV-A71 2A benchmark ``annotated_complexes.csv``
file, selects the benchmark methods/seeds and top-ranked poses, and writes a
standardized table in both CSV and Parquet format.

Supported co-folding input formats are CSV, TSV, and Parquet. Multiple input
files can be supplied with ``--cofolding-files``. By default, outputs are
written to the current working directory; use ``--out-dir`` to choose another
destination.

Example
-------
python util05_Py_prepare_plot_data.py \
    --cofolding-files examples/metrics_boltz2_ost.tsv \
    --annotation-file ../EV-A71_2A_benchmark/structure/processed_outputs/annotated_complexes.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SEEDS = [
    1370180479,
    1449838082,
    1832854922,
    1880307061,
    2012026466,
]

COFOLDING_METHODS = [
    "af3",
    "boltz-1",
    "boltz-2",
    "of3p2",
    "of3p2_ft",
    "protenix",
    "rf3",
]

FINAL_DATA_COLS = [
    "source",
    "method",
    "dock_prot",
    "complex_id",
    "seed",
    "sample",
    "rank",
    "rank_score",
    "rank_by_pair_iptm",
    "pair_iptm",
    "lddt_lp",
    "pocket_bb_rmsd",
    "lig_rmsd",
    "lddt_pli",
    "pb_valid",
    "ligand_smiles",
    "fragment_screen",
    "artefact",
    "pb_valid_groundtruth",
    "filtered",
    "pocket_qcov",
    "sucos_shape",
    "sucos_shape_pocket_qcov",
]


def read_table(path: Path) -> pd.DataFrame:
    """Read a supported co-folding result table."""
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".tsv":
        return pd.read_csv(path, sep="\t")

    raise ValueError(f"Unsupported file type: {path}")


def to_bool(series: pd.Series) -> pd.Series:
    """Convert common boolean representations to bool."""
    if series.dtype == "bool":
        return series.fillna(False)

    return series.fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])


def normalize_seed_column(series: pd.Series) -> pd.Series:
    """Convert seed labels such as seed_123 to integer values."""
    return (
        series.astype(str)
        .str.replace("seed_", "", regex=False)
        .str.extract(r"(\d+)", expand=False)
        .astype("Int64")
    )


def normalize_sample_column(series: pd.Series) -> pd.Series:
    """Extract integer sample identifiers."""
    return series.astype(str).str.extract(r"(\d+)", expand=False).astype("Int64")


def read_annotations(path: Path) -> pd.DataFrame:
    """Read the standard OpenBind complex annotation table."""
    ann = pd.read_csv(path)

    required = {
        "complex_name",
        "smiles",
        "fragment_screen",
        "pb_valid",
        "artefact",
    }
    missing = required - set(ann.columns)
    if missing:
        raise ValueError(
            f"Annotation file is missing required columns: {sorted(missing)}"
        )

    ann = ann.rename(
        columns={
            "complex_name": "complex_id",
            "smiles": "ligand_smiles",
            "pb_valid": "pb_valid_groundtruth",
        }
    )

    ann["fragment_screen"] = to_bool(ann["fragment_screen"])
    ann["pb_valid_groundtruth"] = to_bool(ann["pb_valid_groundtruth"])
    ann["artefact"] = to_bool(ann["artefact"])
    ann["filtered"] = (~ann["pb_valid_groundtruth"]) | ann["artefact"]

    if ann["complex_id"].duplicated().any():
        duplicates = ann.loc[
            ann["complex_id"].duplicated(),
            "complex_id",
        ].unique()
        raise ValueError(f"Duplicate complex_id entries found: {duplicates[:10]}")

    return ann[
        [
            "complex_id",
            "ligand_smiles",
            "fragment_screen",
            "artefact",
            "pb_valid_groundtruth",
            "filtered",
        ]
    ].copy()


def find_rank_column(df: pd.DataFrame) -> str | None:
    """Find the confidence column used to rank co-folding poses."""
    for column in [
        "pair_iptm",
        "ranking_score",
        "confidence_score",
        "iptm",
        "ptm",
    ]:
        if column in df.columns:
            return column

    return None


def read_cofolding_file(path: Path) -> pd.DataFrame:
    """Read one analysed result table and normalize it to the common schema.

    Column aliases used by different co-folding pipelines are normalized,
    invalid poses are removed when ``is_proper`` is available, seed/sample
    identifiers are standardized, and the best available confidence column is
    exposed as ``pair_iptm`` for downstream ranking.
    """
    df = read_table(path)

    df = df.rename(
        columns={
            "target": "complex_id",
            "lddt-pli": "lddt_pli",
            "bb_rmsd_lp": "pocket_bb_rmsd",
            "pocket_qcov-2021-09-30": "pocket_qcov",
            "sucos_shape-2021-09-30": "sucos_shape",
            "sucos_shape_pocket_qcov-2021-09-30": "sucos_shape_pocket_qcov",
        }
    )

    if "is_proper" in df.columns:
        df = df[df["is_proper"]].copy()

    required = {
        "method",
        "complex_id",
        "lig_rmsd",
        "lddt_pli",
        "pb_valid",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

    if "dock_prot" not in df.columns:
        df["dock_prot"] = "cofold"

    if "seed" in df.columns:
        df["seed"] = normalize_seed_column(df["seed"])
    else:
        df["seed"] = 1

    if "sample" in df.columns:
        df["sample"] = normalize_sample_column(df["sample"])
    else:
        df["sample"] = df.groupby(["method", "complex_id", "seed"]).cumcount()

    rank_col = find_rank_column(df)
    if rank_col is None:
        raise ValueError(
            f"No co-folding ranking column found in {path}. "
            "Expected one of: pair_iptm, ranking_score, "
            "confidence_score, iptm, ptm."
        )

    df["pair_iptm"] = pd.to_numeric(df[rank_col], errors="coerce")
    df["lig_rmsd"] = pd.to_numeric(df["lig_rmsd"], errors="coerce")
    df["lddt_pli"] = pd.to_numeric(df["lddt_pli"], errors="coerce")
    df["pb_valid"] = to_bool(df["pb_valid"])
    df["source"] = "cofolding"

    keep_cols = [
        "source",
        "method",
        "dock_prot",
        "complex_id",
        "seed",
        "sample",
        "pair_iptm",
        "lddt_lp",
        "pocket_bb_rmsd",
        "lig_rmsd",
        "lddt_pli",
        "pb_valid",
        "pocket_qcov",
        "sucos_shape",
        "sucos_shape_pocket_qcov",
    ]

    return df[[col for col in keep_cols if col in df.columns]].copy()


def load_cofolding_data(cofolding_files: list[Path]) -> pd.DataFrame:
    """Load and combine explicitly selected co-folding result files.

    Parameters
    ----------
    cofolding_files
        Paths to one or more analysed co-folding result tables. Supported
        formats are CSV, TSV, and Parquet.

    Returns
    -------
    pandas.DataFrame
        Standardized rows from all supplied result files.

    Raises
    ------
    FileNotFoundError
        If any supplied path does not point to an existing file.
    ValueError
        If an input file has an unsupported extension.
    """
    supported_suffixes = {".csv", ".tsv", ".parquet"}
    files = []

    for path in cofolding_files:
        if not path.is_file():
            raise FileNotFoundError(f"Co-folding result file not found: {path}")
        if path.suffix.lower() not in supported_suffixes:
            raise ValueError(f"Unsupported co-folding result file type: {path}")
        files.append(path)

    print(f"Loading {len(files)} co-folding result files")
    return pd.concat(
        [read_cofolding_file(path) for path in files],
        ignore_index=True,
        sort=False,
    )


def add_annotations(
    cofolding_df: pd.DataFrame,
    annotation_file: Path,
) -> pd.DataFrame:
    """Join ground-truth complex annotations onto co-folding pose rows.

    The annotation table is merged many-to-one on ``complex_id`` so every
    predicted pose receives the corresponding ligand and filtering metadata.
    """
    annotations = read_annotations(annotation_file)
    return cofolding_df.merge(
        annotations,
        on="complex_id",
        how="left",
        validate="many_to_one",
    )


def select_final_data(cofolding_df: pd.DataFrame) -> pd.DataFrame:
    """Select benchmark rows and assign per-seed and overall pose ranks.

    Only configured co-folding methods and benchmark seeds are retained. For
    each method/complex/seed combination, the five highest-confidence poses
    are kept. The remaining poses are then ranked across seeds for each
    method/complex pair.
    """
    present_methods = set(cofolding_df["method"].dropna().unique())
    methods = [method for method in COFOLDING_METHODS if method in present_methods]

    df = cofolding_df[
        cofolding_df["method"].isin(methods) & cofolding_df["seed"].isin(SEEDS)
    ].copy()

    if df.empty:
        raise ValueError("No co-folding rows matched the configured methods and seeds.")

    df["rank_score"] = df["pair_iptm"]

    df = df.sort_values(
        ["method", "complex_id", "seed", "pair_iptm", "sample"],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    ).copy()

    df["rank_by_pair_iptm"] = df.groupby(["method", "complex_id", "seed"]).cumcount()

    df = df[df["rank_by_pair_iptm"] < 5].copy()

    df = df.sort_values(
        [
            "method",
            "complex_id",
            "pair_iptm",
            "rank_by_pair_iptm",
            "seed",
            "sample",
        ],
        ascending=[True, True, False, True, True, True],
        kind="mergesort",
    ).copy()

    df["rank"] = df.groupby(["method", "complex_id"]).cumcount()

    for col in FINAL_DATA_COLS:
        if col not in df.columns:
            df[col] = pd.NA

    return (
        df[FINAL_DATA_COLS]
        .sort_values(
            ["method", "complex_id", "rank", "seed", "sample"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


def print_summary(df: pd.DataFrame) -> None:
    """Print basic coverage of the final co-folding table."""
    print(f"\nCo-folding rows: {len(df)}")
    print("\nCo-folding coverage:")
    print(
        df.groupby("method")
        .agg(
            n_rows=("complex_id", "size"),
            n_complexes=("complex_id", "nunique"),
            n_seeds=("seed", "nunique"),
            min_rank=("rank_by_pair_iptm", "min"),
            max_rank=("rank_by_pair_iptm", "max"),
        )
        .reset_index()
        .to_string(index=False)
    )


def write_outputs(df: pd.DataFrame, out_dir: Path) -> None:
    """Write final co-folding CSV and parquet files."""
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "final_cofolding_pose_data.csv"
    parquet_path = out_dir / "final_cofolding_pose_data.parquet"

    df.to_csv(csv_path, index=False)
    df.to_parquet(parquet_path, index=False)

    print("\nWrote final pose-level co-folding data:")
    print(f"CSV:     {csv_path}")
    print(f"Parquet: {parquet_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare final pose-level co-folding results."
    )
    parser.add_argument(
        "--cofolding-files",
        required=True,
        nargs="+",
        type=Path,
        metavar="FILE",
        help=(
            "One or more analysed co-folding result files (.csv, .tsv, or .parquet)."
        ),
    )
    parser.add_argument(
        "--annotation-file",
        required=True,
        type=Path,
        help="Path to annotated_complexes.csv.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("."),
        help=(
            "Directory for final CSV and parquet files (default: current directory)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cofolding_df = load_cofolding_data(args.cofolding_files)
    cofolding_df = add_annotations(
        cofolding_df,
        args.annotation_file,
    )
    final_df = select_final_data(cofolding_df)

    write_outputs(final_df, args.out_dir)
    print_summary(final_df)


if __name__ == "__main__":
    main()

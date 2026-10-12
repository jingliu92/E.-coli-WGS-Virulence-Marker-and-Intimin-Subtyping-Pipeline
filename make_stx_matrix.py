#!/usr/bin/env python3
"""Combine per-assembly StxTyper TSVs into a sample-level matrix.

Expected StxTyper result column: stx_type. Each result file should be named
<SAMPLE_ID>.tsv and be written by the stxtyper command.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

SUFFIXES = (".contigs.fasta", ".fasta", ".fna", ".fa")


def sample_id(path: Path) -> str:
    for suffix in SUFFIXES:
        if path.name.endswith(suffix):
            return path.name[:-len(suffix)]
    return path.stem


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies-dir", default="assembly")
    parser.add_argument("--results-dir", default="stxtyper_results")
    parser.add_argument("--out", default="stx_matrix.tsv")
    args = parser.parse_args()

    assembly_dir, results_dir = Path(args.assemblies_dir), Path(args.results_dir)
    assemblies = sorted({p for pat in ("*.fasta", "*.fa", "*.fna") for p in assembly_dir.glob(pat) if p.is_file()})
    if not assemblies:
        raise SystemExit(f"No assembly FASTA files found in {assembly_dir}")

    rows = []
    for assembly in assemblies:
        sample = sample_id(assembly)
        result_file = results_dir / f"{sample}.tsv"
        if not result_file.exists():
            raise SystemExit(f"Missing StxTyper result for {sample}: {result_file}. Do not interpret a missing result as Stx-negative.")
        try:
            df = pd.read_csv(result_file, sep="\t", comment="#", dtype=str)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()

        if df.empty:
            stx_types = []
        else:
            col_lookup = {c.strip().lower(): c for c in df.columns}
            if "stx_type" not in col_lookup:
                raise SystemExit(f"Expected 'stx_type' column in {result_file}; found columns: {list(df.columns)}")
            stx_types = df[col_lookup["stx_type"]].dropna().astype(str).str.strip().tolist()
            stx_types = [x for x in stx_types if x and x.lower() not in {"nan", "none", "-", "na"}]

        stx1 = sorted({x for x in stx_types if x.lower().startswith("stx1")})
        stx2 = sorted({x for x in stx_types if x.lower().startswith("stx2")})
        rows.append({
            "Sample": sample,
            "stx1": int(bool(stx1)),
            "stx2": int(bool(stx2)),
            "stx1_subtype": ";".join(stx1),
            "stx2_subtype": ";".join(stx2),
        })

    out = Path(args.out)
    pd.DataFrame(rows, columns=["Sample", "stx1", "stx2", "stx1_subtype", "stx2_subtype"]).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out} ({len(rows)} samples)")


if __name__ == "__main__":
    main()

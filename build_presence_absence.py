#!/usr/bin/env python3
"""Build a per-assembly binary marker matrix from filtered BLAST tabular hits."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

FASTA_SUFFIXES = (".contigs.fasta", ".fasta", ".fna", ".fa")


def sample_id(path: Path) -> str:
    name = path.name
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def fasta_ids(path: Path) -> list[str]:
    ids = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    if not ids:
        raise ValueError(f"No FASTA headers found in {path}")
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate marker IDs found in {path}")
    return ids


def read_hit_ids(path: Path) -> set[str]:
    hits = set()
    if not path.exists():
        raise FileNotFoundError(f"Missing BLAST result file: {path}")
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            hits.add(fields[0])  # qseqid = marker ID
    return hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies-dir", default="assembly")
    parser.add_argument("--hits-dir", default="blast_out")
    parser.add_argument("--marker-file", default="markers.fna")
    parser.add_argument("--out", default="blast_presence_absence.tsv")
    args = parser.parse_args()

    assemblies_dir = Path(args.assemblies_dir)
    hits_dir = Path(args.hits_dir)
    assemblies = sorted({p for pat in ("*.fasta", "*.fa", "*.fna") for p in assemblies_dir.glob(pat) if p.is_file()})
    if not assemblies:
        raise SystemExit(f"No FASTA assemblies found in {assemblies_dir}")

    samples = [sample_id(p) for p in assemblies]
    if len(samples) != len(set(samples)):
        raise SystemExit("Assembly filenames produce duplicate sample IDs. Rename files so IDs are unique.")

    markers = fasta_ids(Path(args.marker_file))
    rows = []
    for sample, assembly in zip(samples, assemblies):
        hit_file = hits_dir / f"{sample}.tsv"
        detected = read_hit_ids(hit_file)
        row = {"Sample": sample}
        row.update({marker: int(marker in detected) for marker in markers})
        rows.append(row)

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Sample", *markers], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output} ({len(rows)} samples, {len(markers)} markers)")


if __name__ == "__main__":
    main()

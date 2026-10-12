# E. coli WGS Virulence and Subtyping Pipeline

A Bash/Python workflow for reference-based screening of *Escherichia coli* genome assemblies, Stx typing with StxTyper, optional serotyping with ECTyper, and rule-based isolate grouping.

> **Status:** The repository is assembled from the current project scripts and marker FASTA, with input handling and output integration made consistent. Before biological use or publication, verify the marker sequences and test the workflow against representative known-positive and known-negative assemblies. Software output formats can vary by version.

## What it does

1. Screens each assembly against `markers.fna` with BLASTN.
2. Applies minimum percent-identity and marker query-coverage thresholds.
3. Creates `blast_presence_absence.tsv` for all input assemblies and all marker IDs.
4. Runs StxTyper on each assembly and creates `stx_matrix.tsv`.
5. Optionally runs ECTyper with `--serotype` and creates `serotype.tsv`.
6. Merges results into `combined.tsv`.
7. Creates six rule-based group tables in `Subtype/`.

## Repository files

- `run_ecoli_pipeline.sh` — main workflow.
- `build_presence_absence.py` — creates the BLAST marker matrix.
- `make_stx_matrix.py` — combines per-assembly StxTyper results.
- `markers.fna` — 20 marker/intimin reference sequences supplied for this project.
- `environment.yml` — Conda environment for Python dependencies and BLAST+.
- `docs/workflow.md` — workflow details, assumptions, and limitations.
- `tests/test_marker_headers.py` — basic marker FASTA validation.

## Requirements

- Linux or a compatible Unix-like environment with Bash.
- BLAST+ (`blastn`).
- StxTyper command `stxtyper`, available in the active environment.
- Python 3.10+ and pandas.
- ECTyper only if using optional `--serotype`.

StxTyper and ECTyper are separate tools and may have their own database installation requirements. Install them using their official documentation and ensure their executables are on `PATH` before running this workflow. The script does not automatically activate a named Conda environment; activate the environment that contains your required tools first.

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/ecoli-wgs-virulence-pipeline.git
cd ecoli-wgs-virulence-pipeline
conda env create -f environment.yml
conda activate ecoli-wgs-virulence
```

Install StxTyper in the active environment according to its official instructions. Install ECTyper separately if serotyping is needed. Confirm tool availability:

```bash
blastn -version
stxtyper --help
# Only if using ECTyper:
ectyper --help
```

## Input files

Place one assembled genome per isolate in `assembly/`. Supported extensions are `.fasta`, `.fa`, and `.fna`. Each filename must have a unique sample ID, for example:

```text
assembly/
├── isolate_001.fasta
├── isolate_002.fasta
└── isolate_003.contigs.fasta
```

The sample ID is the filename with the assembly suffix removed. Keep IDs consistent with the names used by downstream tools. Input assemblies and generated outputs are intentionally excluded from Git.

## Run

First check the available options:

```bash
bash run_ecoli_pipeline.sh --help
```

Run the main workflow using default thresholds (90% identity and 90% query coverage):

```bash
bash run_ecoli_pipeline.sh
```

Specify thresholds:

```bash
bash run_ecoli_pipeline.sh --identity 95 --coverage 90
```

Run with ECTyper enabled:

```bash
bash run_ecoli_pipeline.sh --identity 90 --coverage 90 --serotype
```

Use the same sample filenames for all analyses. The workflow requires StxTyper results for every input assembly; if a result file is missing, it stops rather than interpreting a missing result as Stx-negative.

## Outputs

| File/directory | Contents |
|---|---|
| `blast_out/` | Filtered BLAST hit tables, including empty files for assemblies without qualifying hits |
| `blast_presence_absence.tsv` | One row per assembly; one binary column per reference ID |
| `stxtyper_results/` | Per-assembly StxTyper output files |
| `stx_matrix.tsv` | `stx1`, `stx2`, and subtype calls per assembly |
| `serotype.tsv` | ECTyper serotype result when enabled; otherwise header-only table |
| `combined.tsv` | Merged marker, Stx subtype, and optional serotype results |
| `Subtype/eae_stx.tsv` | Stx-positive isolates with `eae` detected |
| `Subtype/stx_subA_no_eae.tsv` | Stx-positive, `subA` positive, `eae` and `aggR` negative |
| `Subtype/stx_aggR_no_eae.tsv` | Stx-positive, `aggR` positive, `eae` and `subA` negative |
| `Subtype/stx_only.tsv` | Stx-positive, `eae`, `subA`, and `aggR` negative |
| `Subtype/no_markers.tsv` | Stx-negative, `eae`, `subA`, and `aggR` negative |
| `Subtype/other_markers.tsv` | Any remaining combination |

The `Subtype/` groups are mutually exclusive due to rule precedence. The first matching group takes precedence; all remaining rows go to `other_markers.tsv`.

## Marker panel

The FASTA contains ten general marker references (`eae`, `stx1_IEH`, `stx2_IEH`, `espK`, `espV`, `espN`, `subA`, `aggR`, `Z2098`, `ureD`) and ten intimin reference sequences (`intimin_gamma`, `intimin_sigma`, `intimin_mu`, `intimin_theta`, `intimin_tau`, `intimin_alpha2`, `intimin_alpha1`, `intimin_beta`, `intimin_rho`, `intimin_epsilon`). The output columns retain the reference IDs as supplied.

The dash in the supplied `espK` sequence was removed because it appeared to be an alignment gap; verify this against the original sequence source before interpreting results. The standalone `A` line within the `Z2098` record is retained as part of the FASTA sequence concatenation; confirm that it is intended. Record reference accession numbers, source publications, and permission/licensing information before public release.

## Interpretation and limitations

- A marker hit is defined by the configured BLAST identity and query-coverage thresholds; it is not proof of expression or phenotype.
- A no-hit result is not necessarily proof of biological absence. Assembly quality, reference choice, sequence divergence, and thresholds can affect detection.
- Intimin reference hits are reported as separate reference-marker columns. This pipeline does **not** automatically assign a definitive intimin subtype based on the best hit.
- `stx1_IEH`/`stx2_IEH` are custom reference-screening markers and are distinct from StxTyper's `stx1`/`stx2` calls.
- The six output groups are operational rules based on `stx1`, `stx2`, `eae`, `subA`, and `aggR`; they are not clinical diagnoses or official pathotype definitions.
- Confirm the StxTyper output schema (`stx_type`) and ECTyper output schema against the installed versions. The pipeline stops when required StxTyper output is missing or lacks the expected column.

## Tests

```bash
python tests/test_marker_headers.py
bash -n run_ecoli_pipeline.sh
bash run_ecoli_pipeline.sh --help
```

A successful syntax/header test does not replace an end-to-end test with representative assemblies and known expected results.

## Citation and license

Cite StxTyper, ECTyper, BLAST+, and the original sources for marker references when using this workflow. Before adding a license, confirm that you have rights to distribute the code and all included reference sequences. The marker sequences may have provenance/licensing requirements separate from the scripts.

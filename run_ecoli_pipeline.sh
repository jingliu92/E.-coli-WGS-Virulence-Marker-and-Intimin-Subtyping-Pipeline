#!/usr/bin/env bash
set -Eeuo pipefail

# E. coli WGS marker screening + StxTyper + optional ECTyper.
IDENTITY=90
COVERAGE=90
ASSEMBLY_DIR="assembly"
MARKER_FILE="markers.fna"
RUN_SEROTYPE=false

usage() {
  cat <<'EOF'
Usage: bash run_ecoli_pipeline.sh [options]

Required input:
  assembly/                  One FASTA assembly per isolate
  markers.fna                Reference marker FASTA

Options:
  --identity PERCENT         Minimum BLAST percent identity (default: 90)
  --coverage PERCENT         Minimum marker query coverage (default: 90)
  --assembly-dir DIR         Assembly directory (default: assembly)
  --marker-file FILE         Marker FASTA (default: markers.fna)
  --serotype                 Run ECTyper (optional; must be installed)
  -h, --help                 Show this help

Required commands for core workflow: blastn, python, stxtyper
Python packages: pandas
EOF
}
while (($#)); do
  case "$1" in
    --identity) [[ $# -ge 2 ]] || { echo "ERROR: --identity needs a value" >&2; exit 2; }; IDENTITY="$2"; shift 2 ;;
    --coverage) [[ $# -ge 2 ]] || { echo "ERROR: --coverage needs a value" >&2; exit 2; }; COVERAGE="$2"; shift 2 ;;
    --assembly-dir) [[ $# -ge 2 ]] || { echo "ERROR: --assembly-dir needs a path" >&2; exit 2; }; ASSEMBLY_DIR="$2"; shift 2 ;;
    --marker-file) [[ $# -ge 2 ]] || { echo "ERROR: --marker-file needs a path" >&2; exit 2; }; MARKER_FILE="$2"; shift 2 ;;
    --serotype) RUN_SEROTYPE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in "$IDENTITY" "$COVERAGE"; do
  awk -v x="$value" 'BEGIN {exit !(x >= 0 && x <= 100)}' || { echo "ERROR: thresholds must be between 0 and 100" >&2; exit 2; }
done
[[ -d "$ASSEMBLY_DIR" ]] || { echo "ERROR: assembly directory not found: $ASSEMBLY_DIR" >&2; exit 1; }
[[ -s "$MARKER_FILE" ]] || { echo "ERROR: marker FASTA not found or empty: $MARKER_FILE" >&2; exit 1; }
command -v blastn >/dev/null || { echo "ERROR: blastn not found. Install BLAST+ first." >&2; exit 1; }
command -v stxtyper >/dev/null || { echo "ERROR: stxtyper not found in PATH. Activate the environment containing StxTyper." >&2; exit 1; }
python -c 'import pandas' 2>/dev/null || { echo "ERROR: pandas is not installed for the active Python." >&2; exit 1; }
if $RUN_SEROTYPE; then command -v ectyper >/dev/null || { echo "ERROR: ectyper not found; install ECTyper or omit --serotype." >&2; exit 1; }; fi

shopt -s nullglob
ASSEMBLIES=("$ASSEMBLY_DIR"/*.fasta "$ASSEMBLY_DIR"/*.fa "$ASSEMBLY_DIR"/*.fna)
shopt -u nullglob
((${#ASSEMBLIES[@]} > 0)) || { echo "ERROR: no .fasta, .fa, or .fna files found in $ASSEMBLY_DIR" >&2; exit 1; }
# Deduplicate filenames and validate unique sample IDs.
declare -A SEEN=()
for f in "${ASSEMBLIES[@]}"; do
  base=$(basename "$f")
  sample="$base"
  for suffix in .contigs.fasta .fasta .fna .fa; do
    if [[ "$sample" == *"$suffix" ]]; then sample="${sample%$suffix}"; break; fi
  done
  [[ -n "$sample" ]] || { echo "ERROR: invalid sample filename: $f" >&2; exit 1; }
  [[ -z "${SEEN[$sample]+x}" ]] || { echo "ERROR: duplicate sample ID '$sample'; rename assembly files." >&2; exit 1; }
  SEEN[$sample]=1
done

mkdir -p blast_out stxtyper_results Subtype
rm -f blast_out/*.tsv

echo "[1/6] BLAST marker screening (identity=${IDENTITY}%, query coverage=${COVERAGE}%)"
for f in "${ASSEMBLIES[@]}"; do
  base=$(basename "$f"); sample="$base"
  for suffix in .contigs.fasta .fasta .fna .fa; do
    if [[ "$sample" == *"$suffix" ]]; then sample="${sample%$suffix}"; break; fi
  done
  raw="blast_out/${sample}.raw.tsv"
  filtered="blast_out/${sample}.tsv"
  blastn -query "$MARKER_FILE" -subject "$f" \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' > "$raw"
  awk -F '\t' -v id="$IDENTITY" -v cov="$COVERAGE" 'BEGIN{OFS="\t"} NF>=11 {qcov=(($7-$6+1)/$5)*100; if($3>=id && qcov>=cov) print $0}' "$raw" > "$filtered"
  rm -f "$raw"
done
python build_presence_absence.py --assemblies-dir "$ASSEMBLY_DIR" --hits-dir blast_out --marker-file "$MARKER_FILE" --out blast_presence_absence.tsv

echo "[2/6] Running StxTyper"
rm -f stxtyper_results/*.tsv
for f in "${ASSEMBLIES[@]}"; do
  base=$(basename "$f"); sample="$base"
  for suffix in .contigs.fasta .fasta .fna .fa; do
    if [[ "$sample" == *"$suffix" ]]; then sample="${sample%$suffix}"; break; fi
  done
  stxtyper -n "$f" --name "$sample" -o "stxtyper_results/${sample}.tsv"
  [[ -f "stxtyper_results/${sample}.tsv" ]] || { echo "ERROR: StxTyper did not create result for $sample" >&2; exit 1; }
done
python make_stx_matrix.py --assemblies-dir "$ASSEMBLY_DIR" --results-dir stxtyper_results --out stx_matrix.tsv

echo "[3/6] Optional ECTyper serotyping"
if $RUN_SEROTYPE; then
  rm -rf Serotype && mkdir -p Serotype
  ectyper -i "${ASSEMBLIES[@]}" -o Serotype
  [[ -s Serotype/output.tsv ]] || { echo "ERROR: expected Serotype/output.tsv from ECTyper" >&2; exit 1; }
  python - <<'PY'
import pandas as pd
from pathlib import Path
p=Path('Serotype/output.tsv')
df=pd.read_csv(p, sep='\t', dtype=str)
if df.empty:
    raise SystemExit('ECTyper output.tsv contains no rows')
cols={c.strip().lower():c for c in df.columns}
sample_col=next((cols[k] for k in ('sample','strain','filename','genome') if k in cols), df.columns[0])
sero_col=next((v for k,v in cols.items() if 'serotype' in k and 'o_type' not in k and 'h_type' not in k), None)
if sero_col is None:
    o_col=next((v for k,v in cols.items() if k in ('o_type','o-type','o antigen')),None)
    h_col=next((v for k,v in cols.items() if k in ('h_type','h-type','h antigen')),None)
    if o_col and h_col:
        df['__serotype']=df[o_col].fillna('') + ':' + df[h_col].fillna('')
        sero_col='__serotype'
    elif len(df.columns)>=8:
        sero_col=df.columns[7]
        print('WARNING: ECTyper has no recognized serotype column; using its eighth column as in the legacy workflow.', flush=True)
    else:
        raise SystemExit(f'Cannot identify a serotype column in ECTyper output. Columns: {list(df.columns)}')
def clean_id(x):
    s=Path(str(x)).name
    for suffix in ('.contigs.fasta','.fasta','.fna','.fa'):
        if s.endswith(suffix): return s[:-len(suffix)]
    return s
out=pd.DataFrame({'Sample':df[sample_col].map(clean_id),'Serotype':df[sero_col].fillna('')})
out.to_csv('serotype.tsv',sep='\t',index=False)
print(f'Wrote serotype.tsv ({len(out)} rows)')
PY
else
  printf 'Sample\tSerotype\n' > serotype.tsv
fi

echo "[4/6] Merging tables"
python - <<'PY'
import pandas as pd
from pathlib import Path
blast=pd.read_csv('blast_presence_absence.tsv',sep='\t',dtype={'Sample':str})
stx=pd.read_csv('stx_matrix.tsv',sep='\t',dtype={'Sample':str})
sero=pd.read_csv('serotype.tsv',sep='\t',dtype={'Sample':str})
for name,df in [('BLAST',blast),('StxTyper',stx),('ECTyper',sero)]:
    if 'Sample' not in df.columns: raise SystemExit(f'{name} table is missing Sample column')
    if df['Sample'].duplicated().any(): raise SystemExit(f'{name} table has duplicate Sample IDs')
combined=blast.merge(stx,on='Sample',how='outer',validate='one_to_one').merge(sero,on='Sample',how='left',validate='one_to_one')
# Keep marker columns (including aggR and all intimin reference hits) and make sure known marker fields exist.
marker_ids=[]
with open('markers.fna') as fh:
    for line in fh:
        if line.startswith('>'): marker_ids.append(line[1:].strip().split()[0])
for col in marker_ids:
    if col not in combined: combined[col]=0
for col in ('stx1','stx2'):
    if col not in combined: combined[col]=0
for col in ('stx1_subtype','stx2_subtype','Serotype'):
    if col not in combined: combined[col]=''
# Fill marker values only; missing Stx calls should be investigated rather than inferred absent.
for col in marker_ids:
    combined[col]=pd.to_numeric(combined[col],errors='coerce').fillna(0).astype(int)
for col in ('stx1','stx2'):
    combined[col]=pd.to_numeric(combined[col],errors='coerce').fillna(0).astype(int)
cols=list(dict.fromkeys(['Sample','stx1_IEH','stx2_IEH','stx1','stx2','stx1_subtype','stx2_subtype',*marker_ids]))
# Add the optional serotype column last.
cols=[c for c in cols if c in combined.columns]+['Serotype']
combined=combined[cols].sort_values('Sample')
combined.to_csv('combined.tsv',sep='\t',index=False)
print(f'Wrote combined.tsv ({len(combined)} samples)')
PY

echo "[5/6] Classifying samples"
python - <<'PYCLASS'
import os
import pandas as pd

df = pd.read_csv('combined.tsv', sep='\t', dtype={'Sample': str})
for col in ('stx1', 'stx2', 'eae', 'subA', 'aggR'):
    if col not in df:
        df[col] = 0
    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
stx_positive = (df['stx1'] == 1) | (df['stx2'] == 1)

candidates = {
    'eae_stx.tsv': stx_positive & (df['eae'] == 1),
    'stx_subA_no_eae.tsv': stx_positive & (df['eae'] == 0) & (df['subA'] == 1) & (df['aggR'] == 0),
    'stx_aggR_no_eae.tsv': stx_positive & (df['eae'] == 0) & (df['subA'] == 0) & (df['aggR'] == 1),
    'stx_only.tsv': stx_positive & (df['eae'] == 0) & (df['subA'] == 0) & (df['aggR'] == 0),
    'no_markers.tsv': (~stx_positive) & (df['eae'] == 0) & (df['subA'] == 0) & (df['aggR'] == 0),
}
assigned = pd.Series(False, index=df.index)
os.makedirs('Subtype', exist_ok=True)
for filename, candidate in candidates.items():
    mask = candidate & ~assigned
    assigned |= mask
    df.loc[mask].to_csv(os.path.join('Subtype', filename), sep='\t', index=False)
    print(f'{filename:<25} {int(mask.sum())} samples')
df.loc[~assigned].to_csv('Subtype/other_markers.tsv', sep='\t', index=False)
print(f"{'other_markers.tsv':<25} {int((~assigned).sum())} samples")
print(f'Total samples: {len(df)}; assigned across groups: {int(assigned.sum() + (~assigned).sum())}')
PYCLASS

echo "[6/6] Complete"
echo "Main outputs: blast_presence_absence.tsv, stx_matrix.tsv, serotype.tsv, combined.tsv"
echo "Classification tables: Subtype/*.tsv"

#!/usr/bin/env bash
# Optional BUSCO / QUAST wrapper. Works on macOS Bash 3.2+. Does not modify references.
set -euo pipefail
if [[ $# -lt 3 ]]; then
  echo 'Usage: bash scripts/03_external_qc.sh config/references.tsv results/02_external_qc /path/to/local_busco_lineage [threads]' >&2; exit 2
fi
config="$1"; out="$2"; lineage="$3"; threads="${4:-4}"
mkdir -p "$out"
command -v busco >/dev/null || { echo 'BUSCO not installed. Core R QC can run without it.' >&2; exit 1; }
[[ -d "$lineage" ]] || { echo "BUSCO lineage directory missing: $lineage" >&2; exit 1; }
{ busco --version || true; command -v quast.py >/dev/null && quast.py --version || true; } > "$out/software_versions.txt" 2>&1
awk -F '\t' 'NR>1 && $1!="" {print $1 "\t" $3 "\t" $5}' "$config" |
while IFS="$(printf '\t')" read -r species genome protein; do
  safe=$(printf '%s' "$species" | tr -c '[:alnum:]_-' '_')
  mkdir -p "$out/$safe"
  busco -i "$genome" -m genome -l "$lineage" -o genome --out_path "$out/$safe" -c "$threads" --offline
  if [[ -n "$protein" && -f "$protein" ]]; then
    busco -i "$protein" -m proteins -l "$lineage" -o proteins --out_path "$out/$safe" -c "$threads" --offline
  fi
  if command -v quast.py >/dev/null 2>&1; then
    quast.py "$genome" -o "$out/$safe/quast" -t "$threads"
  fi
done

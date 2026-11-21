#!/usr/bin/env bash
# Module 1: read-only file inventory and checksums; Bash 3.2+ (macOS-compatible).
set -euo pipefail
if [[ $# -ne 2 ]]; then echo "Usage: bash scripts/01_inventory.sh config/references.tsv results/01_inventory" >&2; exit 2; fi
config="$1"; output="$2"
[[ -f "$config" ]] || { echo "Missing config: $config" >&2; exit 1; }
mkdir -p "$output"
manifest="$output/module1_file_manifest.tsv"
provenance="$output/run_provenance.txt"
printf 'species\taccession\trole\tpath\tsize_bytes\tsha256\n' > "$manifest"
# Portable macOS + Linux hash/size utilities.
if command -v shasum >/dev/null 2>&1; then HASH_TOOL=shasum
elif command -v sha256sum >/dev/null 2>&1; then HASH_TOOL=sha256sum
else echo 'Requires shasum or sha256sum' >&2; exit 1; fi
{
  printf 'run_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf 'config\t%s\n' "$config"
  printf 'host\t%s\n' "$(hostname)"
  printf 'system\t%s\n' "$(uname -s)"
  printf 'hash_method\tSHA256\n'
} > "$provenance"
# Configuration must be tab-separated; tabs and newlines in paths are unsupported.
awk -F '\t' 'NR==1 {if($1!="species" || $2!="accession" || $3!="genome") exit 9; next} NF && NF<5 {exit 10}' "$config" || { echo 'Malformed configuration TSV' >&2; exit 1; }
# Preserve trailing empty fields with per-column extraction in awk, then use tab-delimited records.
# Each role is emitted as its own row, avoiding Bash IFS collapsing empty TSV columns.
awk -F '\t' 'NR>1 && $1!="" {roles[3]="genome";roles[4]="gff3";roles[5]="protein";roles[6]="cds";roles[7]="rna";roles[8]="sequence_report";for(i=3;i<=8;i++) if($i!="") print $1 "\t" $2 "\t" roles[i] "\t" $i}' "$config" |
while IFS="$(printf '\t')" read -r species accession role path; do
  [[ -f "$path" ]] || { echo "Missing $role file for $species: $path" >&2; exit 1; }
  if [[ "$HASH_TOOL" == shasum ]]; then sha=$(shasum -a 256 "$path" | awk '{print $1}')
  else sha=$(sha256sum "$path" | awk '{print $1}'); fi
  size=$(wc -c < "$path" | tr -d ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$species" "$accession" "$role" "$path" "$size" "$sha" >> "$manifest"
done
echo "Created $manifest and $provenance"

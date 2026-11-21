# Mytilus comparative genomics — Modules 1 and 2 (Bash + R)

**No Python required.** This local-first toolkit is intended for macOS/Linux. It does not edit or overwrite reference FASTA, GFF3, or other NCBI source files. It has not been executed against your real NCBI data. The descriptive scripts use Bash and **base R** (no CRAN packages).

## Setup

1. Unzip this package into your working repository, separately from your NCBI reference source directory.
2. Copy `config/references.example.tsv` to `config/references.tsv`.
3. Edit the tab-separated file in a plain text editor or RStudio. For each species, enter the **absolute paths** to `genomic.fna`, `genomic.gff`, `protein.faa`. Other file columns are optional. Some NCBI packages use `genomic.gff.gz` etc.; `.gz` files are supported by the R QC script. If you downloaded all three together, find each species' file in `ncbi_dataset/data/GCF_.../`.
4. Run from the project's top directory:

```bash
bash scripts/01_inventory.sh config/references.tsv results/01_inventory
Rscript scripts/02_genome_qc.R config/references.tsv results/02_qc
Rscript scripts/04_review_tables.R results/02_qc
```

On macOS, RStudio can also run `source('scripts/02_genome_qc.R')` only after setting up `commandArgs`; **recommended execution** is the RStudio Terminal command `Rscript ...` above. Alternatively use `system('Rscript scripts/02_genome_qc.R config/references.tsv results/02_qc')` from the R console.

## Module 1 — provenance / preservation audit

`01_inventory.sh` produces:
- `results/01_inventory/module1_file_manifest.tsv` with species, accession, file role, original path, bytes, and SHA-256 checksum.
- `results/01_inventory/run_provenance.txt` with date, config and platform.

The recorded hash can verify a later transfer by running `shasum -a 256 /new/path/file` and comparing it with the manifest. Local hashes do not independently authenticate NCBI provenance. Do not edit the original source files; add manual acquisition dates, source URLs, BioSample/assembly metadata, and any alternate assemblies to your `docs/decision_log.tsv` or notebook.

**Important:** The inventory hashes the files explicitly listed in the config, not every file in the original unzipped NCBI package. Preserve that package's original README, metadata, and checksum manifest if present.

## Module 2 — descriptive genome and annotation audit

`02_genome_qc.R` produces:
- `module2_assembly_summary.tsv`: FASTA count, total bp, min/max/mean/median length, N50/L50, N90/L90, GC% (among A/C/G/T bases), N% (all sequence positions), invalid sequence characters, zero-length records.
- `module2_annotation_summary.tsv`: gene feature rows, unique gene IDs, missing IDs, explicitly marked protein-coding genes, transcript feature count, CDS and exon features.
- `module2_protein_summary.tsv`: protein FASTA record counts and length summary.
- `module2_gff_feature_counts.tsv`: counts of each GFF3 feature type, for diagnostic review.
- `run_provenance.txt`: R version, configuration and time.

`04_review_tables.R` makes one descriptive four-panel comparison PDF. **A gene feature count is not an inferred number of distinct biological genes across species**; annotation strategies differ. N50 is contiguity, not accuracy. Number of FASTA sequences is **not** the chromosome number. NCBI `sequence_report.jsonl` should be used when you later classify chromosomes vs unplaced/unlocalized scaffolds.

Protein coding depends on the `gene_biotype=` or `gene_type=` attribute in *gene rows*. When neither is present the script reports zero explicitly designated protein-coding gene rows; it does not infer coding status from `CDS` rows. Transcript counts include `mRNA` and `transcript` feature types (a sum of rows, not deduplicated transcript IDs).

## Optional external genome completeness / assembly QC

Install BUSCO and QUAST separately. For reproducibility, obtain one **locally installed** compatible Mollusca BUSCO lineage dataset and run exactly the same release across species. `03_external_qc.sh` requires an existing local lineage folder and runs genome and protein modes plus QUAST (if installed):

```bash
bash scripts/03_external_qc.sh config/references.tsv results/02_external_qc \
  /absolute/path/to/busco_downloads/lineages/mollusca_odb12 4
```

The shown lineage path is **illustrative**: use the actual locally available release consistent with your installed BUSCO version. BUSCO genome mode may use additional gene-prediction software and substantial RAM; run on appropriately provisioned lab hardware if needed. If a BUSCO run fails, do not infer incompleteness—inspect the BUSCO log and dependencies. QUAST outputs are reference-free descriptive metrics (not validated misassembly rates).

## Recommended sequence of decisions

1. Freeze exact assemblies and annotation versions and confirm the selected FASTA/GFF/protein files match each accession.
2. Examine descriptive summary discrepancies; compare against the NCBI assembly data report and sequence report.
3. Compare BUSCO genome mode vs BUSCO proteins mode using the identical dataset/release.
4. Decide whether NCBI annotation is comparable enough for orthology analysis or whether a harmonized reannotation / sensitivity analysis is needed.
5. Record every exclusion, exception and tool version in `docs/decision_log.tsv`.

## Boundaries

No inputs are sent to cloud services. Checksums only read files. Scripts use tab-delimited plain-text config: **do not include tabs or newlines in file paths**, and use absolute paths. Extremely large FASTA/GFF3 files may take time; R reads them in chunks, although sequence lengths and gene IDs remain in RAM. BUSCO is optional and independently installed. Genome BUSCO and proteome BUSCO measure different aspects of completeness; do not combine them into one score.

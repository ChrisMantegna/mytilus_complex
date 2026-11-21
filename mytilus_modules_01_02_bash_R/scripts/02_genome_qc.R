#!/usr/bin/env Rscript
# Module 2: reference genome and annotation descriptive QC; BASE R ONLY.
# Usage: Rscript scripts/02_genome_qc.R config/references.tsv results/02_qc
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript scripts/02_genome_qc.R config/references.tsv results/02_qc')
config_file <- args[1]; out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
refs <- read.delim(config_file, sep = '\t', header = TRUE, quote = '', comment.char = '',
                   check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
required <- c('species', 'accession', 'genome', 'gff3', 'protein')
if (!all(required %in% names(refs))) stop('Config missing columns: ', paste(setdiff(required, names(refs)), collapse = ', '))
if (anyDuplicated(refs$species) || anyDuplicated(refs$accession)) stop('Duplicate species or accession in config')
if (nrow(refs) < 1) stop('No species in config')
if (anyNA(refs[, required])) stop('Required config values cannot be NA')
if (any(!nzchar(refs$genome)) || any(!nzchar(refs$gff3)) || any(!nzchar(refs$protein))) stop('genome, gff3 and protein paths required')
for (p in unique(unlist(refs[, c('genome','gff3','protein')]))) if (!file.exists(p)) stop('File not found: ', p)

open_text <- function(path) {
  if (grepl('\\.gz$', path, ignore.case = TRUE)) gzfile(path, 'rt') else file(path, 'rt')
}
# FASTA parser: streaming; can handle large genome files without loading sequences at once.
# Computes N50/N90 over contig lengths; GC% excludes Ns and other ambiguous bases.
fasta_stats <- function(path, is_protein = FALSE) {
  con <- open_text(path); on.exit(close(con))
  lengths <- numeric(); current <- 0; records <- 0
  bases <- 0; gc <- 0; at <- 0; n_bases <- 0; invalid <- 0
  header_seen <- FALSE
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    for (ln in lines) {
      if (startsWith(ln, '>')) {
        if (header_seen) lengths <- c(lengths, current)
        header_seen <- TRUE; records <- records + 1; current <- 0
      } else if (nzchar(ln)) {
        if (!header_seen) stop('FASTA sequence before header in ', path)
        seq <- toupper(gsub('[[:space:]]+', '', ln))
        s <- nchar(seq, type = 'bytes'); current <- current + s
        if (!is_protein) {
          # Base counting uses regex on each line; handles wrapped sequence lines.
          gc <- gc + nchar(gsub('[^GC]', '', seq))
          at <- at + nchar(gsub('[^AT]', '', seq))
          n_bases <- n_bases + nchar(gsub('[^N]', '', seq))
          invalid <- invalid + nchar(gsub('[ACGTRYSWKMBDHVN.-]', '', seq))
          bases <- bases + s
        }
      }
    }
  }
  if (header_seen) lengths <- c(lengths, current)
  if (!records) stop('No FASTA records in ', path)
  ranked <- sort(lengths, decreasing = TRUE)
  total <- sum(ranked)
  nx <- function(frac) {
    if (total == 0) return(c(N = NA_real_, L = NA_real_))
    idx <- which(cumsum(ranked) >= frac * total)[1]
    c(N = ranked[idx], L = idx)
  }
  n50 <- nx(.5); n90 <- nx(.9)
  data.frame(records = records, total_length = total, min_length = min(lengths),
             max_length = max(lengths), mean_length = mean(lengths),
             median_length = median(lengths), N50 = unname(n50['N']),
             L50 = unname(n50['L']), N90 = unname(n90['N']), L90 = unname(n90['L']),
             gc_percent_acgt = if (!is_protein && (gc + at) > 0) 100*gc/(gc+at) else NA_real_,
             n_percent_all = if (!is_protein && bases > 0) 100*n_bases/bases else NA_real_,
             invalid_characters = if (is_protein) NA_real_ else invalid,
             zero_length_records = sum(lengths == 0), check.names = FALSE)
}
# GFF3 parser: streaming and deliberately conservative about gene type labels.
# Gene ID is read from ID= in gene features. mRNA/other transcript rows counted separately.
gff_stats <- function(path) {
  con <- open_text(path); on.exit(close(con))
  types <- character(); counts <- integer(); gene_ids <- character()
  gene_count <- 0; coding_count <- 0; mrna_count <- 0
  cds_count <- 0; exon_count <- 0; gene_missing_id <- 0
  attr_value <- function(attrs, key) {
    m <- regexpr(paste0('(^|;)', key, '=[^;]*'), attrs, perl = TRUE)
    if (m[1] == -1) return(NA_character_)
    x <- regmatches(attrs, m)
    sub(paste0('^;?', key, '='), '', x)
  }
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, '#')]
    for (ln in lines) {
      fields <- strsplit(ln, '\t', fixed = TRUE)[[1]]
      if (length(fields) < 9) next
      tp <- fields[3]
      pos <- match(tp, types)
      if (is.na(pos)) {types <- c(types, tp); counts <- c(counts, 1L)} else counts[pos] <- counts[pos] + 1L
      if (tp == 'gene') {
        gene_count <- gene_count + 1
        gid <- attr_value(fields[9], 'ID')
        if (is.na(gid) || !nzchar(gid)) gene_missing_id <- gene_missing_id + 1 else gene_ids <- c(gene_ids, gid)
        bio <- attr_value(fields[9], 'gene_biotype')
        if (is.na(bio)) bio <- attr_value(fields[9], 'gene_type')
        if (!is.na(bio) && bio == 'protein_coding') coding_count <- coding_count + 1
      }
      if (tp %in% c('mRNA','transcript')) mrna_count <- mrna_count + 1
      if (tp == 'CDS') cds_count <- cds_count + 1
      if (tp == 'exon') exon_count <- exon_count + 1
    }
  }
  list(summary = data.frame(gene_rows = gene_count, unique_gene_ids = length(unique(gene_ids)),
                             gene_rows_missing_id = gene_missing_id,
                             explicitly_protein_coding_gene_rows = coding_count,
                             transcript_feature_rows = mrna_count,
                             cds_feature_rows = cds_count, exon_feature_rows = exon_count),
       features = data.frame(feature_type = types, count = counts))
}
assemblies <- list(); annotations <- list(); proteins <- list(); features <- list()
for (i in seq_len(nrow(refs))) {
  r <- refs[i, ]; cat('Processing', r$species, '(', r$accession, ')\n')
  ast <- fasta_stats(r$genome)
  gst <- gff_stats(r$gff3)
  pst <- fasta_stats(r$protein, is_protein = TRUE)
  assemblies[[i]] <- cbind(data.frame(species = r$species, accession = r$accession), ast)
  annotations[[i]] <- cbind(data.frame(species = r$species, accession = r$accession), gst$summary)
  proteins[[i]] <- cbind(data.frame(species = r$species, accession = r$accession),
                         pst[, c('records', 'total_length', 'min_length', 'max_length', 'mean_length', 'median_length', 'zero_length_records')])
  features[[i]] <- cbind(data.frame(species = r$species, accession = r$accession), gst$features)
}
write.table(do.call(rbind, assemblies), file.path(out, 'module2_assembly_summary.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
write.table(do.call(rbind, annotations), file.path(out, 'module2_annotation_summary.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
write.table(do.call(rbind, proteins), file.path(out, 'module2_protein_summary.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
write.table(do.call(rbind, features), file.path(out, 'module2_gff_feature_counts.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
writeLines(c(paste0('R version: ', getRversion()),
             paste0('Date UTC: ', format(Sys.time(), tz = 'UTC', usetz = TRUE)),
             paste0('Config: ', normalizePath(config_file))), file.path(out, 'run_provenance.txt'))
cat('Completed Module 2:', normalizePath(out), '\n')

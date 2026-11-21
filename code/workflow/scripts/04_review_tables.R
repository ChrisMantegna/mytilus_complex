#!/usr/bin/env Rscript
# Optional base-R comparison plots; not inferential tests.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: Rscript scripts/04_review_tables.R results/02_qc')
p <- args[1]
a <- read.delim(file.path(p, 'module2_assembly_summary.tsv'), check.names = FALSE)
g <- read.delim(file.path(p, 'module2_annotation_summary.tsv'), check.names = FALSE)
pdf(file.path(p, 'module2_comparison_plots.pdf'), width = 10, height = 7)
par(mfrow = c(2,2), mar = c(7,4,3,1))
barplot(a$total_length/1e9, names.arg = a$species, las = 2, main = 'Assembly length', ylab = 'Gbp')
barplot(a$N50/1e6, names.arg = a$species, las = 2, main = 'N50', ylab = 'Mbp')
barplot(g$gene_rows, names.arg = g$species, las = 2, main = 'Annotated gene features', ylab = 'Count')
barplot(g$explicitly_protein_coding_gene_rows, names.arg = g$species, las = 2,
        main = 'Explicitly protein-coding gene features', ylab = 'Count')
dev.off()
cat('Wrote ', file.path(p, 'module2_comparison_plots.pdf'), '\n', sep = '')

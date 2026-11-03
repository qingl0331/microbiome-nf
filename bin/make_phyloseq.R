#!/usr/bin/env Rscript
# Build a phyloseq object from the amplicon pipeline's exported/ folder.
# Run in your own R (phyloseq, ape installed):
#   Rscript make_phyloseq.R results/exported [sample_metadata.tsv] [out.rds]
suppressPackageStartupMessages({ library(phyloseq); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
ex   <- if (length(args) >= 1) args[1] else "results/exported"
md   <- if (length(args) >= 2 && args[2] != "") args[2] else NA
out  <- if (length(args) >= 3) args[3] else file.path(ex, "phyloseq.rds")

counts <- read.delim(file.path(ex, "feature_table.tsv"), skip = 1, row.names = 1,
                     check.names = FALSE, comment.char = "")
tax <- read.delim(file.path(ex, "taxonomy.tsv"), row.names = 1, check.names = FALSE)
# Name ranks by prefix, not position: current SILVA classifiers add a k__ (kingdom) rank, GTDB / Greengenes2 do not
ranks <- c(d = "Domain", k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
split <- strsplit(tax$Taxon, ";\\s*")
taxm <- t(vapply(split, function(x) {
  out <- setNames(rep(NA_character_, length(ranks)), ranks)
  key <- substr(x, 1, 1)
  val <- sub("^[a-z]__", "", x)
  ok <- grepl("^[a-z]__", x) & key %in% names(ranks) & val != ""
  out[ranks[key[ok]]] <- val[ok]
  out
}, character(length(ranks))))
rownames(taxm) <- rownames(tax)
taxm <- taxm[rownames(counts), colSums(!is.na(taxm)) > 0, drop = FALSE]   # drop ranks this database does not use

ps <- phyloseq(otu_table(as.matrix(counts), taxa_are_rows = TRUE),
               tax_table(taxm),
               phy_tree(read.tree(file.path(ex, "rooted_tree.nwk"))),
               refseq(Biostrings::readDNAStringSet(file.path(ex, "asv_sequences.fasta"))))
if (!is.na(md)) {
  m <- read.delim(md, row.names = 1, check.names = FALSE, comment.char = "")
  m <- m[!grepl("^#q2:types", rownames(m)), , drop = FALSE]
  sample_data(ps) <- sample_data(m)
}
saveRDS(ps, out)
print(ps)

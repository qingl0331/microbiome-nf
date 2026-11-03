# Working with your phyloseq.rds

CBI 16S / ITS amplicon results include `exported/phyloseq.rds`: the whole dataset in one R object, ready for analysis.
It is built by `make_phyloseq.R` in this folder:

```bash
Rscript make_phyloseq.R results/exported [sample_metadata.tsv] [out.rds]   # default out: results/exported/phyloseq.rds
```

Every example below was tested on a phyloseq.rds made by the pipeline (R 4.4, phyloseq from Bioconductor).

## What is inside

phyloseq keeps the parts of a microbiome dataset together. When you remove a sample or merge to genus, every part updates at once.

| Part | What it holds | Get it with |
|---|---|---|
| ASV table | Read counts: one row per ASV (exact sequence variant), one column per sample | `otu_table(ps)` |
| Taxonomy | Name of each ASV, Domain to Genus (and Species when the database has it) | `tax_table(ps)` |
| Tree | Phylogenetic tree of the ASVs, used by UniFrac and Faith's PD | `phy_tree(ps)` |
| Sequences | The DNA sequence of each ASV | `refseq(ps)` |
| Sample data | Your metadata: group, time point, subject… (if you sent it to us) | `sample_data(ps)` |

## 1. Open it in R

You need R with the `phyloseq` package (Bioconductor). One-time install: `BiocManager::install("phyloseq")`.
On CHPC: `module load R/4.4.2`, create `~/R/x86_64-pc-linux-gnu-library/4.4` first, then install.

```r
library(phyloseq)
ps <- readRDS("exported/phyloseq.rds")
ps                       # summary: taxa, samples, ranks, tree

ntaxa(ps); nsamples(ps)
rank_names(ps)           # e.g. Domain Kingdom Phylum Class Order Family Genus
sample_sums(ps)          # reads per sample after all QC
head(tax_table(ps))
```

## 2. Add your sample metadata

Skip this if `ps` already shows `sample_data()`. Otherwise use a tab-separated file whose first column holds the
sample IDs exactly as they appear in the ASV table.

```r
md <- read.delim("metadata.tsv", row.names = 1, check.names = FALSE, comment.char = "")
md <- md[!grepl("^#q2:types", rownames(md)), , drop = FALSE]   # drop the QIIME 2 type row, if any
sample_data(ps) <- sample_data(md)
sample_variables(ps)     # your column names
```

## 3. Remove low-depth samples and rare ASVs

Pick thresholds that suit your study; these are common starting points.

```r
ps <- prune_samples(sample_sums(ps) >= 1000, ps)                    # samples with at least 1,000 reads
ps <- filter_taxa(ps, function(x) sum(x > 0) >= 2, prune = TRUE)   # ASVs seen in at least 2 samples
```

## 4. Look at composition

Convert counts to proportions, merge ASVs to genus, and plot the ten most abundant genera.

```r
ps_rel   <- transform_sample_counts(ps, function(x) x / sum(x))
ps_genus <- tax_glom(ps_rel, "Genus", NArm = FALSE)
top      <- names(sort(taxa_sums(ps_genus), decreasing = TRUE))[1:10]
plot_bar(prune_taxa(top, ps_genus), fill = "Genus")
```

## 5. Alpha diversity: within each sample

How many different ASVs a sample has, and how evenly reads are spread among them.
Replace `group` with a column from your metadata.

```r
plot_richness(ps, x = "group", measures = c("Observed", "Shannon"))
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon"))   # table for statistics
```

## 6. Beta diversity: between samples

How different the communities are. In the PCoA plot, samples with similar communities sit close together.

```r
ord <- ordinate(ps_rel, method = "PCoA", distance = "bray")       # Bray-Curtis
plot_ordination(ps_rel, ord, color = "group")

ord_u <- ordinate(ps, method = "PCoA", distance = "unifrac")      # UniFrac uses the tree
```

## 7. Export tables for Excel or other tools

```r
write.csv(psmelt(ps_genus), "genus_long.csv", row.names = FALSE)    # one row per sample x genus, with metadata
write.csv(as.data.frame(otu_table(ps_genus)), "genus_wide.csv")     # genera x samples
```

## Statistics

Most R microbiome packages take a phyloseq object directly. Choosing the method and model for your design is best
done with Cancer Biostatistics.

| Package | Use |
|---|---|
| DESeq2 | Differential abundance; convert with `phyloseq_to_deseq2()` |
| ANCOM-BC2 | Differential abundance with bias correction; accepts phyloseq |
| MaAsLin2 | Multivariable association; use `otu_table()` and `sample_data()` |
| vegan::adonis2 | PERMANOVA on a distance matrix from `distance()` |

## Good to know

- **Counts are raw.** They are not rarefied or normalised, because each method expects something different.
  Proportions are fine for plots; follow each statistics package's own advice.
- **Mitochondria, chloroplasts, eukaryotes and ASVs without a phylum call are already removed.**
- **`NA` in the taxonomy means "not classified at this rank".** 16S often stops at genus, and Species is absent when
  the reference database has no species level.
- **Rank names depend on the reference database.** SILVA adds a Kingdom rank; GTDB and Greengenes2 do not.
  Check `rank_names(ps)` before using a rank in code.
- **Sample IDs must match exactly** between the ASV table and your metadata, including capitalisation.

Questions about the data processing go to Cancer Bioinformatics through the CBI Jira portal. Questions about study
design and statistics go to Cancer Biostatistics.

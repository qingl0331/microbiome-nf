# CBI microbiome pipelines (v2)

Nextflow (≥ 25.10, classic `params.x` syntax; lint-clean) rewrite of the three microbiome pipelines in `nextflowWorkflow/`:

| Old | New | Status |
|---|---|---|
| `16sPipeline/` (`q2prescreenFilt*.nf` + 3 sbatch scripts, gut/env/EMP variants) | `amplicon.nf` | one workflow, all variants via params |
| `metagenomeHumann3/mg.master.singu.dsl2.hn3.nf` | `shotgun.nf` | DNA samples |
| `microbiomeSingularityPipeline/` (DSL1, HUMAnN2 — no longer runs on Nextflow ≥ 22.12) | `shotgun.nf` | RNA rows in the samplesheet (matched MTX) |
| `q2trainingFeatureClassifiers.sh` | `train_classifier.nf` | one-time per QIIME 2 version |

All tools run from pinned public containers (`conf/containers.config`); no module loads, hard-coded
binary paths or custom Docker images.

## Quick start (CHPC)

**Step-by-step instructions: [`RUNNING.md`](RUNNING.md).** In short:

```bash
module load apptainer nextflow
export NXF_VER=25.10.0                      # tested version; add to ~/.bash_profile
git clone <repo> microbiome-nf && export REPO=$PWD/microbiome-nf
# new site only: edit the EDIT lines in conf/chpc.config, then  sbatch scripts/setup_databases.sh
cd /path/to/project
sbatch $REPO/scripts/run_pipeline.sbatch shotgun.nf  --input samples.csv -output-dir results
sbatch $REPO/scripts/run_pipeline.sbatch amplicon.nf --input samples.csv --classifier <classifier.qza> -output-dir results
```

### 16S / amplicon

```bash
# pass 1 — QC, primer trimming, import, truncation suggestion; stops before DADA2
nextflow run amplicon.nf -profile chpc -work-dir /scratch/general/<fs>/$USER/work_PROJ \
    --input samples.csv --classifier $DB/qiime2_2026.7/classifiers/SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza -output-dir results

less results/trunc/trunc_suggestion.txt          # + results/qiime2/demux/*.qzv on view.qiime2.org

# pass 2 — same command, -resume, add lengths (everything from pass 1 is reused)
nextflow run amplicon.nf ... -resume --trunc_f 250 --trunc_r 210 [--metadata md.tsv --sampling_depth 10000]
```
`--auto_trunc` skips the stop and uses the suggestion directly (logged in `results/trunc/`).

Samplesheet (`assets/example_amplicon_samplesheet.csv`):
```
sample,fastq_1,fastq_2,run
S1,/path/S1_L001_R1_001.fastq.gz,/path/S1_L001_R2_001.fastq.gz,runA
S1,/path/S1_L002_R1_001.fastq.gz,/path/S1_L002_R2_001.fastq.gz,runA     <- same sample = lanes merged
S3,/path/S3_L001_R1_001.fastq.gz,/path/S3_L001_R2_001.fastq.gz,runB     <- DADA2 is run per `run`
```
Or skip the sheet: `--reads '/path/*_R{1,2}_001.fastq.gz'` (sample = text before the first `_`, as before).

Previous variants as params:
- **gut (Greengenes) / environmental (SILVA)**: just a different `--classifier`.
- **EMP V4 protocol**: `--primer_f GTGCCAGCMGCCGCGGTAA --primer_r GGACTACHVGGGTWTCTAAT --amplicon_length 253`
  on per-sample FASTQs (demultiplex first with bcl-convert or `qiime demux emp-paired`).
- **other kits/regions**: `--primer_f/--primer_r/--amplicon_length`.

### Shotgun metagenome / metatranscriptome

```bash
nextflow run shotgun.nf -profile chpc -work-dir /scratch/... --input samples.csv -output-dir results
#   --skip_humann                          taxonomy only (much faster)
#   --host_index /path/mm10.standard   mouse (any bowtie2 index prefix)
#   --phix_index ''                        do not remove PhiX (default: PhiX index in conf/chpc.config)
#   --optical_dedup false                  non-patterned flow cells
```
Samplesheet (`assets/example_shotgun_samplesheet.csv`):
```
sample,fastq_1,fastq_2,type,pair_id
P1_DNA,...,...,DNA,
P1_RNA,...,...,RNA,P1_DNA      <- RNA uses P1_DNA's MetaPhlAn profile for HUMAnN
P3_RNA,...,...,RNA,            <- no matched DNA: MetaPhlAn runs on the RNA itself
```

## What changed and why

### amplicon.nf
| Before | Now | Why |
|---|---|---|
| `cutadapt -u 16 -U 24` (fixed length) | primer trimmed **by sequence**, pairs without both primers discarded, then optional 3' read-through trim | fixed lengths (16/24) did not match the primer lengths in the training script (19/23); sequence trimming also handles staggered primers and short amplicons |
| fastq_screen `--filter` on R1 and R2 **separately** → BBMap repair → Perl barcode re-sync → `exclude-seqs` vs GG 13_8 88_otus | `qiime taxa filter-table/filter-seqs` removing mitochondria, chloroplast, eukaryota and ASVs without a phylum call (`--exclude_taxa`, `--include_taxa`) | filtering mates independently de-pairs reads (the reason repair/re-sync existed); taxonomy-based filtering is the standard approach and needs no host references |
| fastq_screen still available | as a QC **report** only (MultiQC) | contamination % is still useful to clients |
| manual inspection → edit `ptf/ptr` → new sbatch | `SUGGEST_TRUNC` + stop/`-resume` or `--auto_trunc` | suggestion uses median quality, read-length survival (DADA2 drops reads shorter than trunc-len) and the merge-overlap constraint |
| one DADA2 over all samples | DADA2 **per sequencing run**, then exact ASV merge | DADA2 error models are run-specific; runs also denoise in parallel |
| tree and taxonomy computed twice (raw + hits) | once, on the filtered ASVs | halves the slowest steps |
| Greengenes 13_8 / SILVA 132, QIIME 2 2021.11 | QIIME 2 **2026.7**; SILVA / GTDB / Greengenes2 classifiers (`train_classifier.nf`) | GG 13_8 is superseded by Greengenes2; classifiers must match the QIIME 2 scikit-learn version |
| .qza/.qzv only | + `exported/` TSV, BIOM, FASTA, Newick, genus table; `bin/make_phyloseq.R` (lab guide: `bin/README_phyloseq.md`) | direct hand-off to R/phyloseq |
| — | optional `core-metrics-phylogenetic` + alpha rarefaction (`--sampling_depth`) | quick first-look deliverable |

### shotgun.nf
| Before | Now | Why |
|---|---|---|
| kneaddata (Trimmomatic + bowtie2 vs hg38 + PhiX, human hard-coded) | fastp (adapters, Q20 sliding window, **poly-G**, min 50) + bowtie2/samtools host removal with any `--host_index`, then PhiX removal (`--phix_index`, added 2026-10-20) | poly-G tails on NovaSeq were not trimmed; mouse projects need a mouse host; transparent pair handling (both mates must be unmapped) |
| MetaPhlAn run inside HUMAnN (on concatenated reads) | MetaPhlAn run **separately**, paired input, profile passed to HUMAnN | taxonomy-only runs possible; `bowtie2out` kept for cheap re-profiling; RNA can reuse DNA profile |
| `--prescreen-threshold 0.00001` | HUMAnN default (0.01 %) | the old value pulls in nearly every trace pangenome → much slower nucleotide search. Add it back in `conf/modules.config` if there was a specific reason |
| HUMAnN2 + PEAR stitching + SortMeRNA 2.1 (DSL1) | SortMeRNA 4.3.7, paired (`--paired_in`), index built **once** and cached (`storeDir`) | no stitching needed; the old pipeline rebuilt the index for every sample |
| relab tables + STAMP `.spf` | joined RPK, CPM, relab; KO/EC regrouped and named (if `humann_utility_db`); strat/unstrat splits | standard downstream inputs (MaAsLin, MTXmodel) |
| — | `qc/read_tracking.tsv` (pairs per sample after each step) | the first thing clients ask when depth looks low |

### Infrastructure
- **Unknown parameters rejected**: `checkParams()` (`modules/common.nf`) fails on any `--option` not defined at the top of the entry script, so a misspelled parameter (e.g. `--trunclen`) fails immediately instead of being silently ignored. (Classic `params.x = ...` declarations are used instead of the typed `params {}` block, which is unreliable on Nextflow 25.10.0.)
- **Workflow `output {}` block**: fixed, documented results layout; `-output-dir` per project.
- **Retry with more memory** on OOM/timeout kills, capped by `resourceLimits` (no more "requirement exceeds available").
- **Apptainer cache** shared across projects; `/uufs` and `/scratch` bound so config-referenced DBs are visible.
- `tests/run_stub_tests.sh` dry-runs all entry points with `-stub` (no containers or databases needed).

## Pinned versions (`conf/containers.config`)
| Tool | Version | Note |
|---|---|---|
| QIIME 2 | 2026.7 (`quay.io/qiime2/qiime2`) | `taxa filter-table` include/exclude syntax changes in 2026.10 — revisit `FILTER_TAXA` before bumping |
| MetaPhlAn | 4.1.1 + `mpa_vJun23_CHOCOPhlAnSGB_202403` | do **not** upgrade to 4.2 / vJan25 until HUMAnN supports it |
| HUMAnN | 3.9 | stable default. HUMAnN 4 (alpha) needs MetaPhlAn ≤ 4.1.1 + `mpa_vOct22_CHOCOPhlAnSGB_202403`: swap container + `--metaphlan_index` + `-t rel_ab_w_read_stats` in METAPHLAN `ext.args`; output globs already handle 4.x names |
| fastp 1.3.6, bowtie2 2.5.4/samtools 1.21, BBMap 39.18, SortMeRNA 4.3.7, fastq_screen 0.16.0, MultiQC 1.35 | Seqera community images as pinned by nf-core/modules |

## Results layout
```
amplicon:  qc/  trunc/  qiime2/{demux,dada2,artifacts,visualizations,diversity}/  exported/
shotgun:   qc/ (MultiQC, fastp html, read_tracking.tsv)  metaphlan/  metaphlan/per_sample/
           humann_DNA/  humann_RNA/  humann_logs/
pipeline_info/ (timeline, report, trace) in the launch directory
```

## Please verify on first real run
1. **Primers**: defaults are the V3–V4 pair from `q2trainingFeatureClassifiers.sh`. Check the MultiQC cutadapt
   section — if far fewer than ~90 % of pairs pass, the primer sequences or orientation are wrong for that kit.
2. **`--amplicon_length`** (default 460 = upper end of V3–V4 after primer removal) drives the overlap check.
3. **HUMAnN regrouping** (`MERGE_HUMANN`): confirm file names in your `utility_mapping` download
   (`map_ko_uniref90.txt.gz`, `map_level4ec_uniref90.txt.gz`, `map_ko_name.txt.gz`, `map_ec_name.txt.gz`, ...).
4. **Partition sizes** in `conf/chpc.config` (`resourceLimits`) — set to the node type you submit to.
5. Run one previously processed project through both old and new pipelines and compare (ASV counts per
   sample, genus-level composition; MetaPhlAn species correlations).

## Not carried over (easy to add back)
StrainPhlAn (was commented out), STAMP `.spf` exports, SEPP fragment-insertion trees, decontam,
RNA/DNA ratio tables (MTXmodel on the matched `humann_DNA`/`humann_RNA` tables is the current recommendation).

## License
MIT (see `LICENSE`). The tools the pipeline runs (QIIME 2, MetaPhlAn, HUMAnN, etc.) are pulled as containers at run time and keep their own licenses.

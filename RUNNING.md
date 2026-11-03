# Running the CBI microbiome pipelines on CHPC

Step-by-step guide for running `amplicon.nf` (16S/ITS) and `shotgun.nf` (metagenome/metatranscriptome).
See `README.md` for what the pipelines do and `PIPELINE_EXPLANATION.txt` for each step in detail.

Tested on CHPC (hci-kp, Nextflow 25.10.0, Apptainer) in Oct–Nov 2026: all four test tiers pass, and one full
project (20 shotgun samples) was compared with the previous pipeline.

---

## 1. What you need

| | |
|---|---|
| CHPC account | with access to the Slurm account/partition in `conf/chpc.config` (`hci-kp`) |
| Apptainer | `module load apptainer` (all tools run from pinned public containers; no other installs) |
| Nextflow 25.10.0 | `module load nextflow` plus `export NXF_VER=25.10.0` (put the export in `~/.bash_profile`) |
| Databases | already in the core space: `/uufs/chpc.utah.edu/common/home/hcibcore/UNID/microbiomeTestPipe/db` |
| Pipeline code | core copy: `/uufs/chpc.utah.edu/common/home/hcibcore/UNID/microbiomeTestPipe/microbiome-nf`, or `git clone` this repository |

In paths below, `UNID` stands for the uNID of the account that owns the core copy of the pipeline and databases
(the folder under `/uufs/chpc.utah.edu/common/home/hcibcore/`). Below, `REPO` is the pipeline folder:

```bash
export REPO=/uufs/chpc.utah.edu/common/home/hcibcore/UNID/microbiomeTestPipe/microbiome-nf
```

Every run is started with `scripts/run_pipeline.sbatch`, a small Slurm job that runs Nextflow. Nextflow then submits
each pipeline step as its own Slurm job. Do not run Nextflow on a login node.

---

## 2. One-time setup (only for a new installation or a new site)

Skip this section if you use the core copy and its databases.

1. **Edit the `EDIT` lines in `conf/chpc.config`**: `db_root`, `host_index`, `phix_index`, Slurm `queue` and
   `clusterOptions` (account), `resourceLimits` (largest node you submit to) and `apptainer.cacheDir`.
2. **Download the databases** (MetaPhlAn, HUMAnN, SortMeRNA; about 60 GB, one Slurm job):
   `sbatch scripts/setup_databases.sh` (set `DB=` first if `db_root` changed).
3. **16S classifier**: download the SILVA classifier for QIIME 2 2026.7 from
   https://www.arb-silva.de/documentation/classifiers/qiime-2 into `$DB/qiime2_2026.7/classifiers/`.
   For Greengenes2 or GTDB, train one with `train_classifier.nf` (section 6).
4. **Check the installation** (login node, a few minutes; edit `tests/chpc/env.sh` first):
   ```bash
   cd $REPO
   bash tests/chpc/00_preflight.sh --pull     # tools, DBs, host indexes, containers
   bash tests/chpc/01_stub_and_validation.sh  # workflow logic, input validation
   ```
   Both should end with `ALL CHECKS PASSED`. `tests/chpc/README.md` describes tiers 2–4 (real test runs).

---

## 3. 16S / ITS amplicon

### Samplesheet (CSV)

```
sample,fastq_1,fastq_2,run
S1,/path/S1_L001_R1_001.fastq.gz,/path/S1_L001_R2_001.fastq.gz,runA
S1,/path/S1_L002_R1_001.fastq.gz,/path/S1_L002_R2_001.fastq.gz,runA
S3,/path/S3_L001_R1_001.fastq.gz,/path/S3_L001_R2_001.fastq.gz,runB
```

- A sample listed on several rows = several lanes; they are joined.
- `run` = sequencing run. DADA2 learns errors per run, so give each run its own label (optional if there is one run).
- Sample names: letters, numbers, `.`, `_`, `-`; no spaces.
- No samplesheet? `--reads '/path/*_R{1,2}_001.fastq.gz'` (sample name = text before the first `_`).

### Pass 1: QC and truncation suggestion (stops before DADA2)

```bash
mkdir -p /path/to/project && cd /path/to/project
CLS=/uufs/chpc.utah.edu/common/home/hcibcore/UNID/microbiomeTestPipe/db/qiime2_2026.7/classifiers/SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza
sbatch $REPO/scripts/run_pipeline.sbatch amplicon.nf --input samples.csv --classifier $CLS -output-dir results
```

When it finishes, check:
- `results/trunc/trunc_suggestion.txt`: suggested truncation lengths.
- `results/qiime2/demux/*.demux.qzv`: read quality plots (open at https://view.qiime2.org).
- `results/qc/multiqc_report.html`, cutadapt section: about 90% or more of pairs should keep both primers.
  Far fewer means the primer sequences or orientation do not match the kit.

### Pass 2: same command plus the truncation lengths

```bash
sbatch $REPO/scripts/run_pipeline.sbatch amplicon.nf --input samples.csv --classifier $CLS -output-dir results \
    --trunc_f 250 --trunc_r 210 --metadata metadata.tsv --sampling_depth 10000
```

Everything from pass 1 is reused. `--auto_trunc` instead of `--trunc_f/--trunc_r` takes the suggestion directly,
in a single pass.

### Common options

| Option | Default | Use |
|---|---|---|
| `--primer_f`, `--primer_r` | V3–V4 pair from the previous training script | other kits or regions; e.g. EMP V4: `GTGCCAGCMGCCGCGGTAA` / `GGACTACHVGGGTWTCTAAT` |
| `--amplicon_length` | 460 | length after primer removal, for the overlap check (EMP V4: 253) |
| `--metadata` | none | QIIME 2 metadata TSV, first column = sample ID; used in plots and diversity |
| `--sampling_depth` | 0 (off) | > 0 runs core diversity metrics and rarefaction at this depth |
| `--exclude_taxa` | `mitochondria,chloroplast,eukaryota` | removed after classification |
| `--include_taxa` | `p__` | keep only ASVs with a phylum call (`''` to turn off) |

### Results

```
results/qc/                      MultiQC report (cutadapt, fastq_screen)
results/trunc/                   truncation suggestion
results/qiime2/demux|dada2|artifacts|visualizations|diversity/
results/exported/                feature_table.tsv/.biom, taxonomy.tsv, asv_sequences.fasta,
                                 rooted_tree.nwk, genus_table.tsv
```

For R: `Rscript $REPO/bin/make_phyloseq.R results/exported metadata.tsv` writes `phyloseq.rds`; send the lab
`bin/README_phyloseq.md` with it.

---

## 4. Shotgun metagenome / metatranscriptome

### Samplesheet (CSV)

```
sample,fastq_1,fastq_2,type,pair_id
P1_DNA,/path/P1_DNA_R1.fastq.gz,/path/P1_DNA_R2.fastq.gz,DNA,
P1_RNA,/path/P1_RNA_R1.fastq.gz,/path/P1_RNA_R2.fastq.gz,RNA,P1_DNA
P3_RNA,/path/P3_RNA_R1.fastq.gz,/path/P3_RNA_R2.fastq.gz,RNA,
```

- `type`: `DNA` (default) or `RNA`. RNA samples also get rRNA removal (SortMeRNA).
- `pair_id` (RNA only): the matched DNA sample. HUMAnN then uses that sample's MetaPhlAn profile.
- DNA-only projects can use `--reads '/path/*_R{1,2}_001.fastq.gz'` instead.

### Run

```bash
cd /path/to/project
sbatch $REPO/scripts/run_pipeline.sbatch shotgun.nf --input samples.csv -output-dir results
```

### Common options

| Option | Default | Use |
|---|---|---|
| `--skip_humann` | false | taxonomy only (MetaPhlAn); much faster. Re-run without it later; finished steps are reused |
| `--host_index` | human hg38 (`hg38.no_alt`) | bowtie2 index prefix of another host, e.g. mouse `.../Mouse/mm10.standard` |
| `--phix_index` | PhiX | `''` to skip PhiX removal |
| `--optical_dedup` | true | `false` for non-patterned flow cells |
| `--optical_dupedist` | 12000 | 12000 NovaSeq; 2500 HiSeq 3000/4000; 40 NextSeq |

### Results

```
results/qc/                      MultiQC, fastp reports, read_tracking.tsv (read pairs left after each step)
results/metaphlan/               merged and species tables; per_sample/ profiles
results/humann_DNA|humann_RNA/   gene families, KO, EC, pathways (RPK / CPM / relab; stratified + unstratified)
results/humann_logs/
```

Check `qc/read_tracking.tsv` first: it shows how many read pairs were lost to QC, host and rRNA for each sample.

### How long it takes

Measured on project 21319R: 20 oral-swab metagenomes, about 130 M read pairs each, on hci-kp.

- **About 2.5 days** from FASTQs to final tables. Most of that was waiting in the queue: HUMAnN needs 64 GB, so at
  most 6 run at once on hci-kp.
- **About 12 h of computing per sample** (median), mostly HUMAnN. **About 2 h** per sample with `--skip_humann`.

---

## 5. Monitoring, resuming and cleaning up

- **Progress**: `squeue -u $USER`, and the head job log `mbnf.<jobid>.log` in the project folder.
- **Resume**: if a run stops (time limit, node failure, a fixed input), submit the same command again.
- **Reports**: `pipeline_info/` in the project folder holds the Nextflow report, timeline and trace (run time and
  memory of every step).
- **Clean up**: results are real copies, so after checking them delete the work folder
  (`/scratch/general/vast/$USER/mbnf_work/<project folder name>`, or `$WORK` if you set it). Scratch is purged after
  60 days anyway.
- Do not delete `~/.nextflow` while a run is going; Nextflow loads parts of itself from there.

---

## 6. Training a classifier (one-time per QIIME 2 version)

Classifiers must be trained with the same scikit-learn as the QIIME 2 container. Example, Greengenes2 full length:

```bash
DB=/uufs/chpc.utah.edu/common/home/hcibcore/UNID/microbiomeTestPipe/db
sbatch $REPO/scripts/run_pipeline.sbatch train_classifier.nf \
    --ref_seqs 2024.09.backbone.full-length.fna.qza --ref_tax 2024.09.backbone.tax.qza \
    --name gg2_2024.09_full_length -output-dir $DB/qiime2_2026.7/classifiers
```

Add `--primer_f/--primer_r` to train on the amplified region only.

---

## 7. Troubleshooting

| Message | Meaning / fix |
|---|---|
| `Unknown parameter(s): --xyz` | misspelled option; options use underscores (`--trunc_f`) |
| `--classifier is required once truncation lengths are set` | add `--classifier` |
| `No bowtie2 index files found for prefix ...` | wrong `--host_index` / `--phix_index` prefix (no `.1.bt2` suffix) |
| `Provide --input samplesheet.csv or --reads '<glob>'` | no input given |
| Pass 1 keeps far fewer than ~90% of pairs (cutadapt) | primers do not match the kit; set `--primer_f/--primer_r` |
| `qiime2/diversity/diversity_skipped.txt` | fewer than 2 samples reached `--sampling_depth`; lower it |
| A step was killed for memory or time | it is retried automatically with more; if it still fails, raise its label in `conf/base.config` |

Questions: Cancer Bioinformatics Shared Resource (CBI Jira portal).

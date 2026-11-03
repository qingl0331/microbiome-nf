# CHPC test suite

Run everything from the repo root on CHPC (not from the Box copy). Edit `tests/chpc/env.sh` first
(and the `#SBATCH --account/--partition` lines if not `hci-kp`).

| Tier | Command | Where / time | What it proves |
|---|---|---|---|
| 0 | `bash tests/chpc/00_preflight.sh --pull` | login node, minutes (first pull ~30 min) | Nextflow ≥ 25.10, `bin/*` executable, lint, Slurm account, every DB / host index / fastq_screen genome present, every container pulls and contains the tools the pipeline calls |
| 1 | `bash tests/chpc/01_stub_and_validation.sh` | login node, ~10 min | Channel logic for all entry points (`-stub`), also under `-profile chpc`; sample routing (MetaPhlAn skipped for paired RNA, lane merge only for multi-lane samples, skip flags); bad inputs fail fast with clear errors |
| 2 | `sbatch tests/chpc/02_amplicon_mini.sbatch` | Slurm, ~1 h | Real QIIME 2 run on bundled mini data: pass 1 stops before DADA2; pass 2 `-resume` reuses pass 1; DADA2 per run; runs merged; single-sample run; lane merge (S1 = 2× S2); all exports + diversity; `--auto_trunc` + `--reads` glob |
| 3 | `sbatch tests/chpc/make_subsets.sh <real_sheet.csv> 1000000 $TEST_ROOT/subset` then `sbatch tests/chpc/03_shotgun_mini.sbatch` | Slurm, 3–8 h | Real shotgun run on a subsampled real project with production DBs: read tracking consistent, SortMeRNA index cached, MetaPhlAn DB version, `--skip_humann` → full run `-resume` reuses everything, paired RNA uses its DNA profile, HUMAnN KO/EC tables named |
| 4 | full previous project through old and new pipelines | Slurm | Scientific equivalence (see below) |

Each script prints `PASS`/`FAIL` lines and exits non-zero on any failure.
Logs: `$TEST_ROOT/*.log`, Slurm logs `mbnf_*.%j.log` in the repo root, Nextflow `.nextflow.log`.

## Tier 4 — old vs new comparison (manual, one project each)
- **Amplicon**: same raw FASTQs, same classifier DB family. Compare per-sample reads after DADA2
  (`qiime2/dada2/dada2_stats.tsv`), ASV counts, and genus-level relative abundance
  (`exported/genus_table.tsv`): Spearman ρ per sample > 0.9 expected; investigate large shifts in
  Mitochondria/Chloroplast (new filter is taxonomy-based).
- **Shotgun**: compare `read_tracking.tsv` with old kneaddata counts (new should keep slightly
  fewer reads: poly-G + pair-aware host removal), MetaPhlAn species correlation per sample
  (ρ > 0.9), HUMAnN pathway abundances (unstratified, relab) — expect differences from the
  prescreen threshold change (0.00001 → default).

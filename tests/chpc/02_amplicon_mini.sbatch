#!/bin/bash
#SBATCH -J mbnf_t2_amplicon
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH --ntasks=1 --cpus-per-task=2 --mem=8G
#SBATCH --time=12:00:00
#SBATCH -o mbnf_t2_amplicon.%j.log
#
# Tier 2 — REAL amplicon run (containers + Slurm) on the bundled mini data:
#   3 samples / 2 runs, 500 synthetic V3-V4 pairs per file (primers present), S1 = 2 identical lanes, runB = 1 sample.
# This Slurm job is only the Nextflow head; tasks are submitted as their own jobs by -profile chpc.
#   cd microbiome-nf && sbatch tests/chpc/02_amplicon_mini.sbatch          (~30-90 min, mostly queueing)
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source tests/chpc/env.sh
source tests/chpc/lib.sh

W=$TEST_ROOT/work_amp
OUT=$TEST_ROOT/amp
RUN="nextflow run amplicon.nf -profile chpc -work-dir $W --classifier $CLASSIFIER"
# The synthetic reads classify mostly as Eukaryota/Unassigned/no phylum, which the default filter
# (correctly) removes -> empty table. Keep the filter step but only drop mitochondria/chloroplast.
FILT=(--exclude_taxa 'mitochondria,chloroplast' --include_taxa '')
rm -rf "$OUT" "$TEST_ROOT/amp_auto" "$TEST_ROOT/amp_skip"

echo "== pass 1: QC + suggestion, must stop before DADA2"
$RUN --input tests/amplicon_samplesheet.csv -output-dir "$OUT" || fail "pass 1 run"
nonempty "$OUT/trunc/trunc_suggestion.tsv"
nonempty "$OUT/trunc/trunc_suggestion.txt"
for r in runA runB; do qza_ok "$OUT/qiime2/demux/$r.demux.qza"; nonempty "$OUT/qiime2/demux/$r.demux.qzv"; done
nonempty "$OUT/qc/multiqc_report.html"
absent "$OUT/qiime2/dada2"
cat "$OUT/trunc/trunc_suggestion.txt"
read -r TF TR OV ST < <(tail -1 "$OUT/trunc/trunc_suggestion.tsv")
[ "$ST" = "OK" ] && pass "suggestion F=$TF R=$TR overlap=$OV" || fail "suggestion status $ST (F=$TF R=$TR)"
# Primer trimming: >= ~90% of pairs should keep both primers (README 'verify on first run' #1)
mq=$(find -H "$OUT/qc/multiqc_report_data" -name 'multiqc_cutadapt*.txt' | head -1)
if [ -n "$mq" ]; then
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="r_processed")p=i;if($i=="r_written")w=i;if($i=="pairs_processed")p=i;if($i=="pairs_written")w=i}next}
        p&&w{printf "      %s %.1f%% pairs kept\n",$1,100*$w/$p; if($w/$p<0.9)bad=1} END{exit bad}' "$mq" \
        && pass "cutadapt keeps >= 90% of pairs" || fail "cutadapt keeps < 90% of pairs (primer sequence/orientation?)"
fi

echo "== pass 2: -resume with the suggested lengths + diversity"
$RUN --input tests/amplicon_samplesheet.csv -output-dir "$OUT" -resume --trunc_f "$TF" --trunc_r "$TR" --sampling_depth 100 "${FILT[@]}" \
    || fail "pass 2 run"
task_status > "$TEST_ROOT/amp_pass2_tasks.txt"
grep -E '^(CUTADAPT_PRIMERS|QIIME_IMPORT)' "$TEST_ROOT/amp_pass2_tasks.txt" | grep -qv CACHED \
    && fail "pass-1 tasks were re-run instead of CACHED (cache='lenient' / work dir?)" || pass "pass-1 tasks reused (CACHED)"
n_dada=$(grep -c '^DADA2_DENOISE' "$TEST_ROOT/amp_pass2_tasks.txt")
[ "$n_dada" -eq 2 ] && pass "DADA2 ran once per run (2)" || fail "DADA2 tasks = $n_dada, expected 2"

for f in table repseqs taxonomy table.filtered repseqs.filtered rooted_tree; do qza_ok "$OUT/qiime2/artifacts/$f.qza"; done
for f in feature_table.tsv feature_table.biom taxonomy.tsv asv_sequences.fasta rooted_tree.nwk genus_table.tsv; do
    nonempty "$OUT/exported/$f"
done
for f in taxa_barplot taxonomy repseqs; do nonempty "$OUT/qiime2/visualizations/$f.qzv"; done
exists "$OUT/qiime2/diversity/core-metrics"
nonempty "$OUT/qiime2/diversity/alpha_rarefaction.qzv"

# genus table must be collapsed at the g__ rank (SILVA classifiers have an extra k__ rank, so not always level 6)
grep -v '^#' "$OUT/exported/genus_table.tsv" | cut -f1 | grep -q 'g__' && pass "genus_table.tsv collapsed at genus (g__)" \
    || fail "genus_table.tsv rows do not end at g__ (wrong collapse level)"

st=$OUT/qiime2/dada2/dada2_stats.tsv
nonempty "$st"
n_samp=$(grep -vcE '^(sample-id|#q2|id)' "$st")
[ "$n_samp" -eq 3 ] && pass "dada2 stats: 3 samples, 1 header (runs merged)" || fail "dada2 stats has $n_samp sample rows"
# lane merging: S1 = 2 identical lanes -> twice the DADA2 input of S2
awk -F'\t' '$1=="S1"{a=$2}$1=="S2"{b=$2}END{r=a/b; printf "      S1/S2 input = %s/%s\n",a,b; exit !(r>1.8&&r<2.2)}' "$st" \
    && pass "S1 lanes concatenated (~2x S2)" || fail "S1 input is not ~2x S2 (lane merge?)"
n_asv=$(grep -vc '^#' "$OUT/exported/feature_table.tsv")
[ "$n_asv" -gt 0 ] && pass "$n_asv ASVs after taxonomy filter" || fail "0 ASVs after filtering (check FILTER_TAXA / taxonomy.tsv)"
grep -m1 '^#OTU ID' "$OUT/exported/feature_table.tsv" | grep -q 'S3' && pass "single-sample run (runB/S3) present in final table" \
    || fail "S3 (only sample of runB) missing from final table"

echo "== diversity skipped when no 2 samples reach --sampling_depth (run must still succeed)"
$RUN --input tests/amplicon_samplesheet.csv -output-dir "$TEST_ROOT/amp_skip" -resume --trunc_f "$TF" --trunc_r "$TR" \
    --sampling_depth 100000 "${FILT[@]}" || fail "run with --sampling_depth 100000"
nonempty "$TEST_ROOT/amp_skip/qiime2/diversity/diversity_skipped.txt"
absent "$TEST_ROOT/amp_skip/qiime2/diversity/alpha_rarefaction.qzv"
nonempty "$TEST_ROOT/amp_skip/exported/feature_table.tsv"

echo "== auto_trunc + glob input (no samplesheet)"
$RUN --reads 'tests/data/*_R{1,2}_001.fastq.gz' --auto_trunc "${FILT[@]}" -output-dir "$TEST_ROOT/amp_auto" || fail "auto_trunc run"
OUT=$TEST_ROOT/amp_auto
nonempty "$OUT/exported/feature_table.tsv"
grep -m1 '^#OTU ID' "$OUT/exported/feature_table.tsv" | tr '\t' '\n' | tail -n +2 | sort | tr '\n' ' ' | grep -q '^S1 S2 S3 $' \
    && pass "glob mode: samples S1 S2 S3" || fail "glob mode sample names wrong: $(grep -m1 '^#OTU ID' "$OUT/exported/feature_table.tsv")"

summary

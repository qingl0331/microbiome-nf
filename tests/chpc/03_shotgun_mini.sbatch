#!/bin/bash
#SBATCH -J mbnf_t3_shotgun
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH --ntasks=1 --cpus-per-task=2 --mem=8G
#SBATCH --time=72:00:00
#SBATCH -o mbnf_t3_shotgun.%j.log
#
# Tier 3 — REAL shotgun run with the production databases on a subsampled real project.
# The bundled tests/data are random reads (MetaPhlAn/HUMAnN would find nothing), so first:
#   sbatch tests/chpc/make_subsets.sh <real_shotgun_sheet.csv> 1000000 $TEST_ROOT/subset
# with >= 2 DNA samples (one with 2 lanes), 1 RNA with pair_id, 1 RNA without pair_id; set SHOTGUN_SHEET in env.sh.
#   cd microbiome-nf && sbatch tests/chpc/03_shotgun_mini.sbatch           (~3-8 h, HUMAnN dominates)
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source tests/chpc/env.sh
source tests/chpc/lib.sh

[ -s "$SHOTGUN_SHEET" ] || { echo "SHOTGUN_SHEET not found: $SHOTGUN_SHEET (run make_subsets.sh first)"; exit 1; }
W=$TEST_ROOT/work_sg
OUT=$TEST_ROOT/sg
RUN="nextflow run shotgun.nf -profile chpc -work-dir $W --input $SHOTGUN_SHEET -output-dir $OUT"
rm -rf "$OUT"

# expected sets from the samplesheet (header-aware; one entry per sample)
col() { awk -F, -v c="$1" 'NR==1{for(i=1;i<=NF;i++)if($i==c)k=i;next}{print $1","(k?$k:"")}' "$SHOTGUN_SHEET" | sort -u; }
ALL=$(col type | cut -d, -f1)
MPA=$(paste -d, <(col type) <(col pair_id | cut -d, -f2) | awk -F, 'toupper($2)!="RNA" || $3==""{print $1}')
PAIRED_RNA=$(paste -d, <(col type) <(col pair_id | cut -d, -f2) | awk -F, 'toupper($2)=="RNA" && $3!=""{print $1}')
HAS_RNA=$(col type | awk -F, 'toupper($2)=="RNA"' | wc -l)

echo "== run 1: --skip_humann (preprocessing + MetaPhlAn)"
$RUN --skip_humann || fail "skip_humann run"
rt=$OUT/qc/read_tracking.tsv
nonempty "$rt"; nonempty "$OUT/qc/multiqc_report.html"
for s in $ALL; do grep -q "^$s	" "$rt" && pass "read_tracking has $s" || fail "read_tracking lacks $s"; done
# counts must never go up from one step to the next; host_removed should be > 0
awk -F'\t' 'NR>1{p=-1; for(i=3;i<=7;i++){if($i=="NA")continue; if(p>=0 && $i>p){print "      not monotonic: "$0; bad=1} p=$i}
     if($5==0||$5=="NA"){print "      no reads after host removal: "$1; bad=1}} END{exit bad}' "$rt" \
    && pass "read counts monotonic, >0 after host removal" || fail "read_tracking inconsistent"
awk -F'\t' 'NR>1 && $6=="NA"{bad=1} END{exit bad}' "$rt" && pass "PhiX removal step counted for every sample" \
    || fail "phix_removed column empty (is params.phix_index set in conf/chpc.config?)"
column -t "$rt"
if [ "$HAS_RNA" -gt 0 ]; then
    awk -F'\t' 'NR>1 && $2=="RNA" && $7=="NA"{bad=1} END{exit bad}' "$rt" && pass "RNA rows have non_rRNA counts" \
        || fail "RNA sample without SortMeRNA step"
    exists "$(nextflow config -profile chpc -flat . | awk -F" = " '$1=="params.db_cache"{gsub(/\x27/,"",$2);print $2}')/sortmerna_smr_v4.3_default_db/idx"
fi
for s in $MPA; do nonempty "$OUT/metaphlan/per_sample/${s}_metaphlan_profile.tsv"; done
for s in $PAIRED_RNA; do absent "$OUT/metaphlan/per_sample/${s}_metaphlan_profile.tsv"; done
nonempty "$OUT/metaphlan/metaphlan_DNA_merged.tsv"
nonempty "$OUT/metaphlan/metaphlan_DNA_species.tsv"
n_sp=$(grep -c 's__' "$OUT/metaphlan/metaphlan_DNA_species.tsv")
[ "$n_sp" -gt 5 ] && pass "$n_sp species in DNA profiles" || fail "only $n_sp species detected (DB / host removal / depth?)"
[ -z "$(grep -L '^#mpa_vJun23_CHOCOPhlAnSGB_202403' "$OUT"/metaphlan/per_sample/*_metaphlan_profile.tsv)" ] && pass "MetaPhlAn DB = vJun23_202403" \
    || fail "MetaPhlAn profile header does not show vJun23_CHOCOPhlAnSGB_202403"
absent "$OUT/humann_DNA"

echo "== run 2: same command with HUMAnN, -resume (preprocessing + MetaPhlAn must be CACHED)"
$RUN -resume || fail "HUMAnN run"
task_status > "$TEST_ROOT/sg_run2_tasks.txt"
grep -E '^(FASTP|HOST_REMOVE|METAPHLAN|SORTMERNA) ' "$TEST_ROOT/sg_run2_tasks.txt" | grep -qv CACHED \
    && fail "preprocessing/MetaPhlAn re-ran instead of CACHED" || pass "preprocessing + MetaPhlAn reused"
n_h=$(grep -c '^HUMANN ' "$TEST_ROOT/sg_run2_tasks.txt"); n_all=$(echo "$ALL" | wc -w)
[ "$n_h" -eq "$n_all" ] && pass "HUMANN ran for all $n_all samples" || fail "HUMANN tasks = $n_h, expected $n_all"
for s in $PAIRED_RNA; do   # paired RNA must use its DNA partner's MetaPhlAn profile
    d=$(awk -F, -v s="$s" 'NR==1{for(i=1;i<=NF;i++)if($i=="pair_id")k=i;next}$1==s{print $k;exit}' "$SHOTGUN_SHEET")
    wd=$(nextflow log last -f name,workdir | awk -v n="HUMANN ($s)" 'index($0, n) == 1 {print $NF}')
    grep -q -- "--taxonomic-profile ${d}_metaphlan_profile.tsv" "$wd/.command.sh" \
        && pass "$s uses ${d}'s MetaPhlAn profile" || fail "$s HUMAnN not using ${d}'s profile"
done
for t in DNA $( [ "$HAS_RNA" -gt 0 ] && echo RNA ); do
    d=$OUT/humann_$t
    for f in genefamilies_rpk genefamilies_cpm pathabundance pathabundance_relab pathabundance_cpm ko_cpm level4ec_cpm; do nonempty "$d/$f.tsv"; done
    ls "$d"/*_unstratified.tsv >/dev/null 2>&1 && pass "humann_$t strat/unstrat splits" || fail "humann_$t has no *_unstratified.tsv"
    grep -qm1 '^K[0-9]\{5\}: ' "$d/ko_cpm.tsv" && pass "humann_$t KO names attached" \
        || fail "humann_$t ko_cpm.tsv has no named KOs (utility_mapping file names?)"
    awk -F'\t' 'NR>1 && $1=="UNMAPPED"{for(i=2;i<=NF;i++) if($i>900000) bad=1} END{exit bad}' "$d/genefamilies_cpm.tsv" \
        && pass "humann_$t UNMAPPED < 90%" || fail "humann_$t: a sample is >90% UNMAPPED"
done
for s in $ALL; do nonempty "$OUT/humann_logs/$s.humann.log"; done

summary

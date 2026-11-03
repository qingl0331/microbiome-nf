#!/bin/bash
# Tier 1 — workflow logic, no containers / databases / Slurm. Login node, ~5-10 min.
#   a) the repo's own stub suite (tests/run_stub_tests.sh)
#   b) stubs again with -profile chpc (catches clashes between conf/chpc.config params and the entry-script params)
#   c) input validation: bad inputs must fail early with a clear message
set -uo pipefail
source tests/chpc/env.sh
source tests/chpc/lib.sh

W=$TEST_ROOT/work_stub
S=$TEST_ROOT/sheets
mkdir -p "$S"
STUB="-stub -work-dir $W -c tests/chpc/stub_local.config"
SG_DBS="--host_index tests/data/db/host/hg_test --phix_index tests/data/db/phix/phix_test --metaphlan_db tests/data/db/mpa --humann_nucleotide_db tests/data/db/choco
        --humann_protein_db tests/data/db/uniref --humann_utility_db tests/data/db/util
        --sortmerna_ref tests/data/db/rrna.fasta --db_cache $TEST_ROOT/stub_cache"

echo "== a) repo stub suite"
# run_stub_tests.sh writes to /tmp; send it to scratch instead
sed "s#/tmp/#$TEST_ROOT/stub_#g" tests/run_stub_tests.sh > "$TEST_ROOT/run_stub_tests.sh"
bash "$TEST_ROOT/run_stub_tests.sh" >"$TEST_ROOT/stub_suite.log" 2>&1 \
    && pass "tests/run_stub_tests.sh" || fail "tests/run_stub_tests.sh (see $TEST_ROOT/stub_suite.log)"

echo "== b) stubs under -profile chpc"
OUT=$TEST_ROOT/stub_chpc_amp
nextflow run amplicon.nf -profile chpc $STUB --input tests/amplicon_samplesheet.csv --classifier tests/data/dummy_classifier.qza \
    --trunc_f 250 --trunc_r 210 --sampling_depth 100 -output-dir "$OUT" >"$TEST_ROOT/stub_chpc_amp.log" 2>&1 \
    && pass "amplicon stub, chpc profile" || fail "amplicon stub, chpc profile (see $TEST_ROOT/stub_chpc_amp.log)"
nonempty "$OUT/trunc/trunc_suggestion.tsv"
exists "$OUT/exported/feature_table.tsv"
exists "$OUT/qiime2/diversity/alpha_rarefaction.qzv"
for r in runA runB; do exists "$OUT/qiime2/demux/$r.demux.qza"; done

OUT=$TEST_ROOT/stub_chpc_sg
nextflow run shotgun.nf -profile chpc $STUB --input tests/shotgun_samplesheet.csv $SG_DBS -output-dir "$OUT" \
    >"$TEST_ROOT/stub_chpc_sg.log" 2>&1 && pass "shotgun stub, chpc profile" || fail "shotgun stub, chpc profile (see $TEST_ROOT/stub_chpc_sg.log)"
# MetaPhlAn must run on DNA samples + the unpaired RNA (P3_RNA) only; P1_RNA borrows P1_DNA's profile
for s in P1_DNA P2_DNA P3_RNA; do exists "$OUT/metaphlan/per_sample/${s}_metaphlan_profile.tsv"; done
absent "$OUT/metaphlan/per_sample/P1_RNA_metaphlan_profile.tsv"
exists "$OUT/humann_DNA"; exists "$OUT/humann_RNA"
n_humann=$(task_status | grep -c '^HUMANN (')
[ "$n_humann" -eq 4 ] && pass "HUMANN ran for all 4 samples (incl. paired RNA)" || fail "HUMANN tasks = $n_humann, expected 4"
n_lanes=$(task_status | grep -c '^MERGE_LANES')
[ "$n_lanes" -eq 1 ] && pass "MERGE_LANES only for the 2-lane sample" || fail "MERGE_LANES tasks = $n_lanes, expected 1"

nextflow run shotgun.nf -profile chpc $STUB --input tests/shotgun_samplesheet.csv $SG_DBS --skip_humann --optical_dedup false \
    -output-dir "$TEST_ROOT/stub_chpc_sg2" >"$TEST_ROOT/stub_chpc_sg2.log" 2>&1 \
    && pass "shotgun stub --skip_humann --optical_dedup false" || fail "shotgun --skip_humann (see $TEST_ROOT/stub_chpc_sg2.log)"
task_status | grep -qE '^(HUMANN|CLUMPIFY)' && fail "HUMANN/CLUMPIFY ran despite skip flags" || pass "skip flags respected"
ls "$W"/*/*/P1_DNA_phix.bowtie2.log >/dev/null 2>&1 && pass "PhiX pass ran in HOST_REMOVE" || fail "no P1_DNA_phix.bowtie2.log in the HOST_REMOVE work dirs"

# PhiX pass is optional: an empty --phix_index must still run (host removal only)
nextflow run shotgun.nf -profile chpc $STUB --input tests/shotgun_samplesheet.csv $SG_DBS --phix_index '' --skip_humann \
    -output-dir "$TEST_ROOT/stub_chpc_sg3" >"$TEST_ROOT/stub_chpc_sg3.log" 2>&1 \
    && pass "shotgun stub without PhiX index" || fail "shotgun stub without PhiX index (see $TEST_ROOT/stub_chpc_sg3.log)"

echo "== c) input validation (each must FAIL fast)"
A0="nextflow run amplicon.nf -profile local $STUB -output-dir $TEST_ROOT/stub_bad"
A="$A0 --classifier tests/data/dummy_classifier.qza"
G="nextflow run shotgun.nf -profile local $STUB $SG_DBS -output-dir $TEST_ROOT/stub_bad"

expect_fail "misspelled param --trunclen"        "Unknown parameter"  $A --input tests/amplicon_samplesheet.csv --trunclen 250
expect_fail "no input"                          "Provide --input"  $A
expect_fail "trunc lengths without classifier"  "--classifier is required" $A0 --input tests/amplicon_samplesheet.csv --trunc_f 250 --trunc_r 210
expect_fail "classifier file missing"           "" $A0 --input tests/amplicon_samplesheet.csv --trunc_f 250 --trunc_r 210 --classifier /nope.qza

printf 'sample,fastq_1,fastq_2\nS 1,tests/data/S1_L001_R1_001.fastq.gz,tests/data/S1_L001_R2_001.fastq.gz\n' > "$S/bad_name.csv"
expect_fail "sample name with space"            "must start with" $A --input "$S/bad_name.csv"
printf 'sample,fastq_1,fastq_2\nS1,tests/data/MISSING_R1.fastq.gz,tests/data/S1_L001_R2_001.fastq.gz\n' > "$S/missing.csv"
expect_fail "FASTQ does not exist"              "" $A --input "$S/missing.csv"
printf 'sample,fastq_1\nS1,tests/data/S1_L001_R1_001.fastq.gz\n' > "$S/nocol.csv"
expect_fail "missing fastq_2 column"            "missing column" $A --input "$S/nocol.csv"
printf 'sample,fastq_1,fastq_2,run\nS1,tests/data/S1_L001_R1_001.fastq.gz,tests/data/S1_L001_R2_001.fastq.gz,runA\nS1,tests/data/S1_L002_R1_001.fastq.gz,tests/data/S1_L002_R2_001.fastq.gz,runB\n' > "$S/conflict.csv"
expect_fail "same sample in two runs"           "different" $A --input "$S/conflict.csv"

printf 'sample,fastq_1,fastq_2,type,pair_id\nX_RNA,tests/data/S3_L001_R1_001.fastq.gz,tests/data/S3_L001_R2_001.fastq.gz,RNA,NOT_THERE\n' > "$S/bad_pair.csv"
expect_fail "RNA pair_id not a DNA sample"      "is not a DNA sample" $G --input "$S/bad_pair.csv"
printf 'sample,fastq_1,fastq_2,type,pair_id\nX,tests/data/S3_L001_R1_001.fastq.gz,tests/data/S3_L001_R2_001.fastq.gz,WGS,\n' > "$S/bad_type.csv"
expect_fail "type not DNA/RNA"                  "type must be DNA or RNA" $G --input "$S/bad_type.csv"
expect_fail "host index prefix without files"   "No bowtie2 index files" $G --input tests/shotgun_samplesheet.csv --host_index /nope/hg38
expect_fail "PhiX index prefix without files"   "phix_index" $G --input tests/shotgun_samplesheet.csv --phix_index /nope/phix
expect_fail "RNA rows but no --sortmerna_ref"   "sortmerna_ref is required" \
    nextflow run shotgun.nf -profile local $STUB --input tests/shotgun_samplesheet.csv --host_index tests/data/db/host/hg_test \
    --metaphlan_db tests/data/db/mpa --skip_humann -output-dir $TEST_ROOT/stub_bad

summary

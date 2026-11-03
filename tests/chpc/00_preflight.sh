#!/bin/bash
# Tier 0 — environment preflight. Login node, ~2 min (+ container pulls the first time).
#   bash tests/chpc/00_preflight.sh            checks only
#   bash tests/chpc/00_preflight.sh --pull     also pulls/converts every pinned container into apptainer.cacheDir
set -uo pipefail
source tests/chpc/env.sh
source tests/chpc/lib.sh

echo "== tools"
ver=$(nextflow -version 2>/dev/null | grep -oE 'version [0-9]+\.[0-9]+\.[0-9]+' | awk '{print $2}')
if [ -n "$ver" ] && [ "$(printf '%s\n25.10.0\n' "$ver" | sort -V | head -1)" = "25.10.0" ]; then
    pass "nextflow $ver (>= 25.10)"
else
    fail "nextflow >= 25.10 required (found '${ver:-none}'): curl -s https://get.nextflow.io | bash"
fi
command -v apptainer >/dev/null && pass "apptainer $(apptainer --version)" || fail "apptainer not on PATH"
command -v sbatch >/dev/null && pass "slurm" || fail "sbatch not on PATH"

echo "== repo"
for f in bin/*.py bin/*.R; do [ -x "$f" ] && pass "executable $f" || fail "$f not executable (Box/zip copies drop +x): chmod +x bin/*"; done
nextflow lint . >"$TEST_ROOT/lint.log" 2>&1 && pass "nextflow lint" || fail "nextflow lint (see $TEST_ROOT/lint.log)"
nextflow config -profile chpc . >"$TEST_ROOT/config_chpc.txt" 2>&1 && pass "config resolves with -profile chpc" \
    || fail "nextflow config -profile chpc (see $TEST_ROOT/config_chpc.txt)"

echo "== slurm account/partition"
sacctmgr -nP show assoc user="$USER" format=account,partition 2>/dev/null | grep -q "$ACCOUNT" \
    && pass "$USER has account $ACCOUNT" || fail "no association for account $ACCOUNT (check conf/chpc.config queue/clusterOptions)"
sinfo -h -p "$PARTITION" -o '%c %m %l' 2>/dev/null | head -1 | while read -r c m l; do
    echo "      $PARTITION node: ${c} cpus, $((m / 1024)) GB, walltime $l  -> compare with resourceLimits in conf/chpc.config"
done

echo "== databases"
exists "$DB/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202403/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl"
exists "$DB/humann3.9/chocophlan"
exists "$DB/humann3.9/uniref"
exists "$DB/humann3.9/utility_mapping"
for m in map_ko_uniref90.txt.gz map_level4ec_uniref90.txt.gz map_ko_name.txt.gz map_ec_name.txt.gz; do
    [ -n "$(find "$DB/humann3.9/utility_mapping" -name "$m" 2>/dev/null | head -1)" ] && pass "utility_mapping has $m" \
        || fail "utility_mapping lacks $m (MERGE_HUMANN expects it directly in that dir; KO/EC tables will be skipped)"
done
exists "$DB/sortmerna/smr_v4.3_default_db.fasta"
exists "$CLASSIFIER"
host=$(nextflow config -profile chpc -flat . 2>/dev/null | awk -F" = " '$1=="params.host_index"{gsub(/\x27/,"",$2);print $2}')
ls "${host}"*.bt2* >/dev/null 2>&1 && pass "host index $host" || fail "no bowtie2 files for host_index '$host'"
phix=$(nextflow config -profile chpc -flat . 2>/dev/null | awk -F" = " '$1=="params.phix_index"{gsub(/\x27/,"",$2);print $2}')
ls "${phix}"*.bt2* >/dev/null 2>&1 && pass "PhiX index $phix" || fail "no bowtie2 files for phix_index '$phix'"
while read -r _ name path; do
    ls "${path}"*.bt2* >/dev/null 2>&1 || fail "fastq_screen DB $name missing: $path"
done < <(grep -P '^DATABASE\t' assets/fastq_screen.conf)

if [ "${1:-}" = "--pull" ]; then
    echo "== containers (pull + probe)"
    cache=$(nextflow config -profile chpc -flat . | awk -F" = " '$1=="apptainer.cacheDir"{gsub(/\x27/,"",$2);print $2}')
    export APPTAINER_CACHEDIR=${cache}/.build
    grep -oE "container = '[^']+'" conf/containers.config | cut -d"'" -f2 | sort -u | while read -r img; do
        sif="$cache/$(echo "$img" | sed 's#[/:]#-#g').img"   # same naming Nextflow uses
        [ -s "$sif" ] || apptainer pull --name "$(basename "$sif")" --dir "$cache" "docker://$img" >/dev/null 2>&1
        [ -s "$sif" ] && pass "image $img" || { fail "pull failed: $img"; continue; }
        case "$img" in
            *qiime2*)        probe="qiime info && cutadapt --version && biom --version && python3 -c 'import sklearn'" ;;
            *fastq-screen*)  probe="fastq_screen --version && bowtie2 --version" ;;   # fastq_screen needs bowtie2 IN this image
            *bowtie2*)       probe="bowtie2 --version && samtools --version && pigz --version" ;;
            *bbmap*)         probe="clumpify.sh --version; pigz --version" ;;
            *fastp*)         probe="fastp --version" ;;
            *sortmerna*)     probe="sortmerna --version" ;;
            *metaphlan*)     probe="metaphlan --version && merge_metaphlan_tables.py -h" ;;
            *humann*)        probe="humann --version && humann_join_tables -h && humann_regroup_table -h" ;;
            *multiqc*)       probe="multiqc --version && python3 --version" ;;   # also READ_TRACKING
            *)               probe="true" ;;
        esac
        # pixi-based Wave images: same PATH as containerOptions in conf/containers.config
        opts=(); apptainer exec "$sif" test -e /shell-hook.sh 2>/dev/null \
            && opts=(--env "PATH=/opt/wave/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin")
        apptainer exec "${opts[@]}" "$sif" bash -c "$probe" >"$TEST_ROOT/probe.log" 2>&1 \
            && pass "  tools present: ${probe//&&/,}" || fail "  probe failed in $img: $(tail -1 "$TEST_ROOT/probe.log")"
    done
fi
summary

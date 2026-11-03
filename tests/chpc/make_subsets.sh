#!/bin/bash
#SBATCH -J mbnf_subset
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH --ntasks=1 --cpus-per-task=4 --mem=16G
#SBATCH --time=4:00:00
#SBATCH -o mbnf_subset.%j.log
#
# Subsample a REAL samplesheet (amplicon or shotgun) to N read pairs per FASTQ pair, with the
# pipeline's own BBMap container. All other columns (run / type / pair_id) are kept as-is.
#   sbatch tests/chpc/make_subsets.sh <samplesheet.csv> <pairs_per_file> <out_dir>
# e.g. shotgun: pick 3 DNA + 1 RNA with pair_id + 1 RNA without, 1,000,000 pairs
#      amplicon: pick ~8 samples from 2 runs, 20,000 pairs
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source tests/chpc/env.sh

SHEET=${1:?samplesheet}; N=${2:?pairs per file}; DEST=${3:?output dir}
mkdir -p "$DEST"
cache=$(nextflow config -profile chpc -flat . | awk -F" = " '$1=="apptainer.cacheDir"{gsub(/\x27/,"",$2);print $2}')
IMG=$(grep -oE "community.wave.seqera.io/library/bbmap_pigz:[^']+" conf/containers.config)
SIF=$cache/$(echo "$IMG" | sed 's#[/:]#-#g').img
[ -s "$SIF" ] || apptainer pull --dir "$cache" --name "$(basename "$SIF")" "docker://$IMG"

out_sheet=$DEST/$(basename "${SHEET%.csv}")_subset.csv
head -1 "$SHEET" | tr -d '\r' > "$out_sheet"
i=0
tail -n +2 "$SHEET" | tr -d '\r' | while IFS=, read -r sample r1 r2 rest; do
    [ -z "$sample" ] && continue
    i=$((i + 1))
    o1=$DEST/${sample}_row${i}_R1.fastq.gz; o2=$DEST/${sample}_row${i}_R2.fastq.gz
    if [ ! -s "$o2" ]; then
        # per-row seed: a sample repeated on several rows (mock lanes) gets different reads
        apptainer exec -B /uufs/chpc.utah.edu -B /scratch "$SIF" reformat.sh -Xmx12g \
            in1="$r1" in2="$r2" out1="$o1" out2="$o2" samplereadstarget="$N" sampleseed=$((41 + i)) ow=t 2>>"$DEST/reformat.log"
    fi
    echo "$sample,$o1,$o2${rest:+,$rest}" >> "$out_sheet"
done
echo "Wrote $out_sheet"; cat "$out_sheet"

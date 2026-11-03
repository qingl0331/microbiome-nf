#!/bin/bash
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH -J mb_db_setup
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH -o mb_db_setup.%j.log
#
# One-time download of versioned databases for the CBI microbiome pipelines.
# Paths match conf/chpc.config. Uses the SAME containers as the pipeline, so the database
# versions are guaranteed to match the tool versions. Safe to re-run (skips finished parts).
#
# Approximate sizes: MetaPhlAn vJun23 ~20 GB, ChocoPhlAn ~16 GB, UniRef90 DIAMOND ~20 GB,
# utility mapping ~2.5 GB, SortMeRNA ~0.1 GB.
set -euo pipefail

DB=${DB:-/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db}   # = params.db_root
module load apptainer 2>/dev/null || module load singularity
RUN="apptainer exec -B /uufs/chpc.utah.edu -B /scratch"

MPA_IMG=docker://quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0
HUMANN_IMG=docker://quay.io/biocontainers/humann:3.9--py312hdfd78af_0
MPA_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403

# ---- MetaPhlAn 4.1.1 marker DB (the version HUMAnN 3.9 expects) ----
if [ ! -e "$DB/metaphlan/$MPA_INDEX/$MPA_INDEX.pkl" ]; then
    mkdir -p "$DB/metaphlan/$MPA_INDEX"
    $RUN $MPA_IMG metaphlan --install --index $MPA_INDEX --bowtie2db "$DB/metaphlan/$MPA_INDEX" --nproc 8
fi

# ---- HUMAnN 3.9 databases ----
H=$DB/humann3.9
mkdir -p "$H"
[ -d "$H/chocophlan" ]      || $RUN $HUMANN_IMG humann_databases --download chocophlan full "$H" --update-config no
[ -d "$H/uniref" ]          || $RUN $HUMANN_IMG humann_databases --download uniref uniref90_diamond "$H" --update-config no
[ -d "$H/utility_mapping" ] || $RUN $HUMANN_IMG humann_databases --download utility_mapping full "$H" --update-config no

# ---- SortMeRNA rRNA reference (metatranscriptomes) ----
if [ ! -e "$DB/sortmerna/smr_v4.3_default_db.fasta" ]; then
    mkdir -p "$DB/sortmerna" && cd "$DB/sortmerna"
    wget -q https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
    tar -xzf database.tar.gz && rm database.tar.gz
    cd - >/dev/null
fi

# ---- QIIME 2 2026.7 classifiers ----
# SILVA: download the pre-trained classifier matching QIIME 2 >= 2026.4 from
#        https://www.arb-silva.de/documentation/classifiers/qiime-2
# Greengenes2 / anything else: train once with
#        nextflow run train_classifier.nf -profile chpc --ref_seqs ... --ref_tax ... --name ...
#        -output-dir $DB/qiime2_2026.7/classifiers
mkdir -p "$DB/qiime2_2026.7/classifiers"

echo "Done. Record versions:"
$RUN $MPA_IMG metaphlan --version
$RUN $HUMANN_IMG humann --version

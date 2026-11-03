#!/bin/bash
#SBATCH -J mbnf
#SBATCH --account=hci-kp
#SBATCH --partition=hci-kp
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH -o mbnf.%j.log
#
# Nextflow head job for the CBI microbiome pipelines. Only this small controller runs here; every
# pipeline step is submitted as its own Slurm job by -profile chpc. Submit from the project folder:
#
#   cd /path/to/project
#   sbatch $REPO/scripts/run_pipeline.sbatch amplicon.nf --input samples.csv --classifier X.qza -output-dir results
#   sbatch $REPO/scripts/run_pipeline.sbatch shotgun.nf  --input samples.csv -output-dir results
#
# Submitting the same command again resumes the run (-resume): finished steps are reused.
# Overrides: REPO (pipeline folder), WORK (Nextflow work dir), NXF_VER (Nextflow version).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

REPO=${REPO:-/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/microbiome-nf}
WORK=${WORK:-/scratch/general/vast/$USER/mbnf_work/$(basename "$PWD")}
if [ $# -lt 1 ]; then
    echo "usage: sbatch run_pipeline.sbatch <amplicon.nf|shotgun.nf|train_classifier.nf> [pipeline options]" >&2
    exit 1
fi
ENTRY=$1
shift

module load apptainer
command -v nextflow >/dev/null || module load nextflow
export NXF_VER=${NXF_VER:-25.10.0}          # the version the pipeline was tested with
export NXF_OPTS='-Xms1g -Xmx4g' NXF_ANSI_LOG=false

echo "START $(date '+%F %T')  $ENTRY $*"
echo "      pipeline: $REPO   work dir: $WORK"
nextflow run "$REPO/$ENTRY" -profile chpc -work-dir "$WORK" -resume "$@"
echo "END   $(date '+%F %T')"

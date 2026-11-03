# Shared settings for the CHPC test suite. EDIT the marked lines once, then:
#   source tests/chpc/env.sh        (every script below does this itself)

# ---------- EDIT ----------
export ACCOUNT=hci-kp
export PARTITION=hci-kp
export DB=/uufs/chpc.utah.edu/common/home/hcibcore/u0762203/microbiomeTestPipe/db   # = params.db_root
export CLASSIFIER=$DB/qiime2_2026.7/classifiers/SILVA_144_SSURef_NR99_uniform_classifier_V3V4-341f-806r.qza                      # real QIIME 2 2026.7 classifier
export TEST_ROOT=/scratch/general/vast/$USER/mbnf_test                                # scratch file system
# Shotgun mini test: samplesheet written by make_subsets.sh (subset of a real, previously processed project)
export SHOTGUN_SHEET=$TEST_ROOT/subset/shotgun_subset.csv
# --------------------------

export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false
mkdir -p "$TEST_ROOT"

module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
command -v nextflow >/dev/null || module load nextflow 2>/dev/null || true

# All scripts must run from the repo root (relative paths in tests/*.csv)
if [ ! -f amplicon.nf ]; then
    echo "Run from the microbiome-nf repo root (cd there before running / sbatch-ing)" >&2
    exit 1
fi

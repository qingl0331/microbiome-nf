#!/bin/bash
# Dry-run every entry point with -stub (no containers/databases needed). Run from repo root.
set -euo pipefail
C="-profile local -stub -ansi-log false"
nextflow lint .
nextflow run amplicon.nf $C --input tests/amplicon_samplesheet.csv --classifier tests/data/dummy_classifier.qza -output-dir /tmp/t_amp1
nextflow run amplicon.nf $C --input tests/amplicon_samplesheet.csv --classifier tests/data/dummy_classifier.qza --trunc_f 250 --trunc_r 210 --sampling_depth 1000 -resume -output-dir /tmp/t_amp1
nextflow run amplicon.nf $C --reads 'tests/data/*_R{1,2}_001.fastq.gz' --classifier tests/data/dummy_classifier.qza --auto_trunc -output-dir /tmp/t_amp2
nextflow run shotgun.nf $C --input tests/shotgun_samplesheet.csv --host_index tests/data/db/host/hg_test --phix_index tests/data/db/phix/phix_test \
  --metaphlan_db tests/data/db/mpa --humann_nucleotide_db tests/data/db/choco --humann_protein_db tests/data/db/uniref \
  --humann_utility_db tests/data/db/util --sortmerna_ref tests/data/db/rrna.fasta --db_cache /tmp/t_cache -output-dir /tmp/t_sg
echo "All stub tests passed"

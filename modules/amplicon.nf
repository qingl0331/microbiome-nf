/*
 * Amplicon (16S/ITS) processes — QIIME 2 2026.7 container for everything.
 */

// Reverse complement including IUPAC degenerate bases
def revcomp(String s) {
    def comp = [A: 'T', C: 'G', G: 'C', T: 'A', U: 'A', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
                B: 'V', V: 'B', D: 'H', H: 'D', N: 'N', I: 'I']
    return s.toUpperCase().reverse().collect { c -> comp[c] ?: 'N' }.join('')
}

// Keep QIIME 2 / matplotlib / numba caches inside the task dir (container home is read-only)
def q2env() {
    return '''
    export TMPDIR="$PWD/tmp" XDG_CONFIG_HOME="$PWD/.config" MPLCONFIGDIR="$PWD/.mpl" NUMBA_CACHE_DIR="$PWD/.numba" HOME="$PWD"
    mkdir -p "$TMPDIR" "$XDG_CONFIG_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
    '''.stripIndent()
}

/*
 * Primer removal by sequence (not fixed-length -u/-U), so staggered/phased primers,
 * degenerate positions and read-through on short amplicons are all handled.
 *   pass 1: 5' primers required on BOTH mates (pairs lacking either are discarded)
 *   pass 2: optional 3' read-through trimming (reverse-complemented opposite primer)
 * Multiple lanes/runs of the same library are concatenated first.
 */
process CUTADAPT_PRIMERS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(r1s, stageAs: 'in_R1/*'), path(r2s, stageAs: 'in_R2/*')
    val primer_f
    val primer_r
    val min_len

    output:
    tuple val(meta), path("${meta.id}_R1.trim.fastq.gz"), path("${meta.id}_R2.trim.fastq.gz"), emit: reads
    path "${meta.id}.cutadapt.log", emit: log

    script:
    def rc_f = revcomp(primer_f)
    def rc_r = revcomp(primer_r)
    """
    # named after the sample: MultiQC takes the cutadapt sample name from these file names
    cat in_R1/* > ${meta.id}_R1.fastq.gz
    cat in_R2/* > ${meta.id}_R2.fastq.gz

    cutadapt -j ${task.cpus} \\
        -g ${primer_f} -G ${primer_r} \\
        --discard-untrimmed --pair-filter any \\
        -e 0.1 --overlap 10 \\
        -o p1_R1.fastq.gz -p p1_R2.fastq.gz \\
        ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz > ${meta.id}.cutadapt.log

    cutadapt -j ${task.cpus} \\
        -a ${rc_r} -A ${rc_f} \\
        -e 0.1 --overlap 10 -m ${min_len} \\
        -o ${meta.id}_R1.trim.fastq.gz -p ${meta.id}_R2.trim.fastq.gz \\
        p1_R1.fastq.gz p1_R2.fastq.gz > readthrough.log

    rm -f ${meta.id}_R?.fastq.gz p1_R?.fastq.gz
    """

    stub:
    """
    echo "" | gzip > ${meta.id}_R1.trim.fastq.gz
    echo "" | gzip > ${meta.id}_R2.trim.fastq.gz
    touch ${meta.id}.cutadapt.log
    """
}

/*
 * Suggest DADA2 truncation lengths from the trimmed reads: per-position median quality,
 * the fraction of reads long enough to survive truncation, and the merge-overlap constraint.
 */
process SUGGEST_TRUNC {
    label 'process_single'

    input:
    path r1s, stageAs: 'r1/*'
    path r2s, stageAs: 'r2/*'
    val amplicon_length
    val min_overlap
    val q_threshold

    output:
    path "trunc_suggestion.tsv", emit: tsv
    path "trunc_suggestion.txt", emit: report

    script:
    """
    suggest_trunc.py \\
        --r1 r1/* --r2 r2/* \\
        --amplicon-length ${amplicon_length} \\
        --min-overlap ${min_overlap} \\
        --q ${q_threshold} \\
        --out-tsv trunc_suggestion.tsv > trunc_suggestion.txt
    cat trunc_suggestion.txt
    """

    stub:
    """
    printf 'trunc_f\\ttrunc_r\\texpected_overlap\\tstatus\\n270\\t220\\t30\\tOK\\n' > trunc_suggestion.tsv
    touch trunc_suggestion.txt
    """
}

/*
 * Build a manifest and import one sequencing run; demux summary for manual inspection.
 */
process QIIME_IMPORT {
    tag "${run}"
    label 'process_low'

    input:
    tuple val(run), val(ids), path(r1s), path(r2s)

    output:
    tuple val(run), path("${run}.demux.qza"), emit: qza
    path "${run}.demux.qzv", emit: qzv

    script:
    def rows = [ids, r1s, r2s].transpose().collect { id, a, b -> "${id}\t\$PWD/${a}\t\$PWD/${b}" }.join('\n')
    """
    ${q2env()}
    printf 'sample-id\\tforward-absolute-filepath\\treverse-absolute-filepath\\n' > manifest.tsv
    cat >> manifest.tsv <<EOF
${rows}
EOF

    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-format PairedEndFastqManifestPhred33V2 \\
        --input-path manifest.tsv \\
        --output-path ${run}.demux.qza

    qiime demux summarize --i-data ${run}.demux.qza --o-visualization ${run}.demux.qzv
    """

    stub:
    """
    touch ${run}.demux.qza ${run}.demux.qzv
    """
}

/*
 * DADA2 per sequencing run (error models are run-specific), threads = allocated CPUs.
 */
process DADA2_DENOISE {
    tag "${run}"
    label 'process_high'

    input:
    tuple val(run), path(demux), val(trunc_f), val(trunc_r)
    val max_ee_f
    val max_ee_r

    output:
    tuple val(run), path("${run}.table.qza"), path("${run}.repseqs.qza"), emit: results
    path "${run}.dada2_stats.tsv", emit: stats
    path "${run}.dada2_stats.qzv", emit: qzv

    script:
    def args = task.ext.args ?: ''
    """
    ${q2env()}
    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs ${demux} \\
        --p-trunc-len-f ${trunc_f} \\
        --p-trunc-len-r ${trunc_r} \\
        --p-max-ee-f ${max_ee_f} \\
        --p-max-ee-r ${max_ee_r} \\
        --p-n-threads ${task.cpus} \\
        ${args} \\
        --output-dir dada2

    mv dada2/table.qza ${run}.table.qza
    mv dada2/representative_sequences.qza ${run}.repseqs.qza
    mv dada2/denoising_stats.qza ${run}.dada2_stats.qza

    qiime metadata tabulate --m-input-file ${run}.dada2_stats.qza --o-visualization ${run}.dada2_stats.qzv
    qiime tools export --input-path ${run}.dada2_stats.qza --output-path stats_export
    mv stats_export/stats.tsv ${run}.dada2_stats.tsv
    """

    stub:
    """
    touch ${run}.table.qza ${run}.repseqs.qza ${run}.dada2_stats.qzv
    printf 'sample-id\\tinput\\n' > ${run}.dada2_stats.tsv
    """
}

/*
 * Merge per-run feature tables / ASV sequences (ASVs are exact sequences, so merging is exact).
 */
process MERGE_RUNS {
    label 'process_low'

    input:
    path tables, stageAs: 'tables/*'
    path seqs, stageAs: 'seqs/*'
    path stats, stageAs: 'stats/*'

    output:
    path "table.qza", emit: table
    path "repseqs.qza", emit: repseqs
    path "dada2_stats.tsv", emit: stats

    script:
    """
    ${q2env()}
    n=\$(ls tables/*.qza | wc -l)
    if [ "\$n" -eq 1 ]; then
        cp tables/*.qza table.qza
        cp seqs/*.qza repseqs.qza
    else
        qiime feature-table merge --i-tables tables/*.qza --o-merged-table table.qza
        qiime feature-table merge-seqs --i-data seqs/*.qza --o-merged-data repseqs.qza
    fi

    # one stats table across runs (keep a single header + q2:types line)
    first=1
    for f in stats/*.tsv; do
        if [ \$first -eq 1 ]; then cat "\$f"; first=0; else grep -v -e '^sample-id' -e '^#q2:types' -e '^id' "\$f" || true; fi
    done > dada2_stats.tsv
    """

    stub:
    """
    touch table.qza repseqs.qza dada2_stats.tsv
    """
}

process CLASSIFY {
    label 'process_medium'
    label 'process_high_memory'

    input:
    path repseqs
    path classifier

    output:
    path "taxonomy.qza", emit: taxonomy

    script:
    """
    ${q2env()}
    qiime feature-classifier classify-sklearn \\
        --i-classifier ${classifier} \\
        --i-reads ${repseqs} \\
        --p-n-jobs ${task.cpus} \\
        --o-classification taxonomy.qza
    """

    stub:
    """
    touch taxonomy.qza
    """
}

/*
 * Replaces fastq_screen --filter + BBMap repair + barcode re-sync + exclude-seqs(88_otus):
 * drop ASVs classified as mitochondria/chloroplast/eukaryota and ASVs without a phylum call.
 */
process FILTER_TAXA {
    label 'process_low'

    input:
    path table
    path repseqs
    path taxonomy
    val exclude_terms
    val include_terms

    output:
    path "table.filtered.qza", emit: table
    path "repseqs.filtered.qza", emit: repseqs

    script:
    def exc = exclude_terms ? "--p-exclude '${exclude_terms}'" : ''
    def inc = include_terms ? "--p-include '${include_terms}'" : ''
    """
    ${q2env()}
    qiime taxa filter-table --i-table ${table} --i-taxonomy ${taxonomy} ${exc} ${inc} \\
        --o-filtered-table table.filtered.qza
    qiime taxa filter-seqs --i-sequences ${repseqs} --i-taxonomy ${taxonomy} ${exc} ${inc} \\
        --o-filtered-sequences repseqs.filtered.qza
    """

    stub:
    """
    touch table.filtered.qza repseqs.filtered.qza
    """
}

process PHYLOGENY {
    label 'process_medium'

    input:
    path repseqs

    output:
    path "rooted_tree.qza", emit: rooted
    path "unrooted_tree.qza", emit: unrooted

    script:
    """
    ${q2env()}
    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences ${repseqs} \\
        --p-n-threads ${task.cpus} \\
        --output-dir phylogeny
    mv phylogeny/rooted_tree.qza rooted_tree.qza
    mv phylogeny/tree.qza unrooted_tree.qza
    """

    stub:
    """
    touch rooted_tree.qza unrooted_tree.qza
    """
}

/*
 * Visualizations for the client / for QC (all from the final filtered table).
 */
process QIIME_VISUALIZE {
    label 'process_low'

    input:
    path table_raw
    path table
    path repseqs
    path taxonomy
    path metadata

    output:
    path "*.qzv", emit: qzv

    script:
    """
    ${q2env()}
    qiime feature-table summarize --i-table ${table_raw} --output-dir summary_unfiltered
    qiime feature-table summarize --i-table ${table} --output-dir summary_filtered
    for d in summary_unfiltered summary_filtered; do
        for v in \$d/*.qzv; do mv "\$v" "\${d}_\$(basename \$v)"; done
    done
    qiime feature-table tabulate-seqs --i-data ${repseqs} --i-taxonomy ${taxonomy} --o-visualization repseqs.qzv \\
      || qiime feature-table tabulate-seqs --i-data ${repseqs} --o-visualization repseqs.qzv
    qiime metadata tabulate --m-input-file ${taxonomy} --o-visualization taxonomy.qzv
    qiime taxa barplot --i-table ${table} --i-taxonomy ${taxonomy} --m-metadata-file ${metadata} \\
        --o-visualization taxa_barplot.qzv
    """

    stub:
    """
    touch taxa_barplot.qzv taxonomy.qzv repseqs.qzv summary_filtered_summary.qzv
    """
}

/*
 * Plain-text exports for R (phyloseq/mia) and for clients who do not use QIIME 2.
 */
process EXPORT_TABLES {
    label 'process_low'

    input:
    path table
    path repseqs
    path taxonomy
    path tree

    output:
    path "exported/*", emit: files

    script:
    """
    ${q2env()}
    mkdir -p exported
    qiime tools export --input-path ${table} --output-path ex_table
    biom convert -i ex_table/feature-table.biom -o exported/feature_table.tsv --to-tsv
    mv ex_table/feature-table.biom exported/feature_table.biom
    qiime tools export --input-path ${taxonomy} --output-path ex_tax
    mv ex_tax/taxonomy.tsv exported/taxonomy.tsv
    qiime tools export --input-path ${repseqs} --output-path ex_seqs
    mv ex_seqs/dna-sequences.fasta exported/asv_sequences.fasta
    qiime tools export --input-path ${tree} --output-path ex_tree
    mv ex_tree/tree.nwk exported/rooted_tree.nwk

    # genus-level counts, a common client request. The genus level is wherever "g__" sits in the taxonomy
    # strings: 6 for GTDB / Greengenes2, 7 for the current SILVA classifiers (which add a k__ kingdom rank).
    glevel=\$(awk -F'\t' 'NR > 1 { n = split(\$2, r, /; */); for (i = 1; i <= n; i++) if (r[i] ~ /^g__/) { print i; exit } }' exported/taxonomy.tsv)
    qiime taxa collapse --i-table ${table} --i-taxonomy ${taxonomy} --p-level \${glevel:-6} --o-collapsed-table genus.qza
    qiime tools export --input-path genus.qza --output-path ex_genus
    biom convert -i ex_genus/feature-table.biom -o exported/genus_table.tsv --to-tsv
    """

    stub:
    """
    mkdir -p exported
    touch exported/feature_table.tsv exported/taxonomy.tsv exported/asv_sequences.fasta exported/rooted_tree.nwk
    """
}

process DIVERSITY {
    label 'process_medium'

    input:
    path table
    path tree
    path metadata
    val depth

    output:
    path "core-metrics", optional: true, emit: core
    path "alpha_rarefaction.qzv", optional: true, emit: rarefaction
    path "diversity_skipped.txt", optional: true, emit: skipped

    script:
    """
    ${q2env()}
    # core-metrics drops samples below the sampling depth and beta diversity needs >= 2 left:
    # skip (with a note) instead of failing the run and cancelling the exports.
    n=\$(python3 -c "import qiime2, biom; t = qiime2.Artifact.load('${table}').view(biom.Table); print(int((t.sum('sample') >= ${depth}).sum()))")
    if [ "\$n" -lt 2 ]; then
        echo "Diversity skipped: \$n sample(s) with >= ${depth} reads after filtering (need 2). Lower --sampling_depth?" > diversity_skipped.txt
        exit 0
    fi
    qiime diversity core-metrics-phylogenetic \\
        --i-phylogeny ${tree} --i-table ${table} \\
        --p-sampling-depth ${depth} \\
        --m-metadata-file ${metadata} \\
        --p-n-jobs-or-threads ${task.cpus} \\
        --output-dir core-metrics
    qiime diversity alpha-rarefaction \\
        --i-table ${table} --i-phylogeny ${tree} \\
        --p-max-depth ${depth} \\
        --m-metadata-file ${metadata} \\
        --o-visualization alpha_rarefaction.qzv
    """

    stub:
    """
    mkdir core-metrics
    touch alpha_rarefaction.qzv
    """
}

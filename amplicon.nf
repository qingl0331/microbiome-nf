#!/usr/bin/env nextflow
/*
 * CBI 16S / amplicon pipeline (replaces 16sPipeline/: q2prescreenFilt*.nf + q2pipe16s*.sh)
 *
 *   primer trim (cutadapt, by sequence) -> [fastq_screen QC] -> QIIME 2 import per run
 *   -> truncation suggestion -> DADA2 per run -> merge runs -> classify-sklearn
 *   -> taxonomy-based off-target filter -> MAFFT/FastTree -> qzv + TSV/BIOM/Newick exports
 *   -> [core diversity metrics]
 *
 * Two-pass use (same as before, but one command and fully cached):
 *   1) nextflow run amplicon.nf -profile chpc --input samples.csv --classifier X.qza
 *        -> stops after QC; read results/trunc/trunc_suggestion.txt and results/qiime2/*.demux.qzv
 *   2) same command + -resume --trunc_f 250 --trunc_r 210     (or use --auto_trunc in step 1)
 */

include { checkParams; optStr; readSamplesheet; readGlob; FASTQ_SCREEN; MULTIQC } from './modules/common.nf'
include {
    CUTADAPT_PRIMERS; SUGGEST_TRUNC; QIIME_IMPORT; DADA2_DENOISE; MERGE_RUNS; CLASSIFY;
    FILTER_TAXA; PHYLOGENY; QIIME_VISUALIZE; EXPORT_TABLES; DIVERSITY
} from './modules/amplicon.nf'

// ---- input (one of) ----
params.input            = null    // CSV: sample,fastq_1,fastq_2[,run]  (repeat a sample on several rows to merge lanes)
params.reads            = null    // or a glob, e.g. '/path/*_R{1,2}_001.fastq.gz'; sample = text before first '_'
params.metadata         = null    // QIIME 2 sample metadata TSV (first column = sample ID); optional

// ---- primers (defaults = V3-V4 set from q2trainingFeatureClassifiers.sh; confirm against the kit) ----
params.primer_f         = 'TGCCTACGGGNBGCASCAG'
params.primer_r         = 'GCGACTACNVGGGTATCTAATCC'
params.min_len          = 50

// ---- DADA2 ----
params.trunc_f          = 0       // 0 = not set
params.trunc_r          = 0
params.auto_trunc       = false   // use the suggested lengths without stopping
params.amplicon_length  = 460     // post-primer amplicon length used for the overlap check (V3-V4 upper range)
params.min_overlap      = 20
params.trunc_q          = 25
params.max_ee_f         = 2.0
params.max_ee_r         = 2.0

// ---- taxonomy ----
params.classifier       = null    // naive Bayes .qza trained with THIS QIIME 2 version (see train_classifier.nf)
params.exclude_taxa     = 'mitochondria,chloroplast,eukaryota'
params.include_taxa     = 'p__'   // require at least a phylum-level call ('' to disable)

// ---- optional ----
params.sampling_depth   = 0       // >0 runs core-metrics-phylogenetic + alpha rarefaction
params.fastq_screen_conf = null   // contamination screen report (set in conf/chpc.config)

workflow {
    main:
    // ------------------------------------------------------------------ inputs
    checkParams(['input', 'reads', 'metadata', 'primer_f', 'primer_r', 'min_len', 'trunc_f', 'trunc_r',
                 'auto_trunc', 'amplicon_length', 'min_overlap', 'trunc_q', 'max_ee_f', 'max_ee_r',
                 'classifier', 'exclude_taxa', 'include_taxa', 'sampling_depth', 'fastq_screen_conf'])
    if (!params.input && !params.reads) {
        error("Provide --input samplesheet.csv or --reads '<glob>'")
    }
    def samples = params.input
        ? readSamplesheet(file(params.input, checkIfExists: true), ['run'])
        : readGlob(params.reads as String, [run: ''])
    samples = samples.collect { meta, r1, r2 -> [meta + [run: meta.run ?: 'run1'], r1, r2] }

    def trunc_given = params.trunc_f > 0 && params.trunc_r > 0
    def will_denoise = trunc_given || params.auto_trunc
    if (will_denoise && !params.classifier) {
        error("--classifier is required once truncation lengths are set (--trunc_f/--trunc_r or --auto_trunc)")
    }
    log.info("""
        CBI amplicon pipeline
        samples   : ${samples.size()} (${samples.collect { s -> s[0].run }.unique().size()} run(s))
        primers   : F=${params.primer_f} R=${params.primer_r}
        truncation: ${trunc_given ? "${params.trunc_f}/${params.trunc_r}" : (params.auto_trunc ? 'auto' : 'not set -> stop after QC + suggestion')}
        classifier: ${params.classifier ?: '-'}
        """.stripIndent())

    ch_samples = channel.fromList(samples)

    // ------------------------------------------------------------------ primer trimming + QC
    CUTADAPT_PRIMERS(ch_samples, params.primer_f, params.primer_r, params.min_len)
    ch_trimmed = CUTADAPT_PRIMERS.out.reads

    ch_screen = channel.empty()
    if (params.fastq_screen_conf) {
        FASTQ_SCREEN(ch_trimmed.map { meta, r1, _r2 -> [meta, r1] }, file(params.fastq_screen_conf, checkIfExists: true))
        ch_screen = FASTQ_SCREEN.out.txt
    }

    SUGGEST_TRUNC(
        ch_trimmed.map { _meta, r1, _r2 -> r1 }.collect(),
        ch_trimmed.map { _meta, _r1, r2 -> r2 }.collect(),
        params.amplicon_length, params.min_overlap, params.trunc_q
    )

    // ------------------------------------------------------------------ import per run
    ch_by_run = ch_trimmed
        .map { meta, r1, r2 -> [meta.run, meta.id, r1, r2] }
        .groupTuple()
    QIIME_IMPORT(ch_by_run)

    MULTIQC(
        CUTADAPT_PRIMERS.out.log.mix(ch_screen).collect(),
        'CBI amplicon QC'
    )

    // ------------------------------------------------------------------ truncation decision
    def ch_trunc
    if (trunc_given) {
        ch_trunc = channel.value([params.trunc_f, params.trunc_r])
    }
    else if (params.auto_trunc) {
        ch_trunc = SUGGEST_TRUNC.out.tsv
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                if (row.status != 'OK') {
                    log.warn("Auto truncation status=${row.status} — review trunc_suggestion.txt")
                }
                [row.trunc_f as Integer, row.trunc_r as Integer]
            }
            .first()
    }
    else {
        ch_trunc = channel.empty()
    }

    // ------------------------------------------------------------------ DADA2 -> taxonomy -> tree
    DADA2_DENOISE(
        QIIME_IMPORT.out.qza.combine(ch_trunc),   // -> [run, demux.qza, trunc_f, trunc_r]
        params.max_ee_f, params.max_ee_r
    )
    MERGE_RUNS(
        DADA2_DENOISE.out.results.map { _run, t, _s -> t }.collect(),
        DADA2_DENOISE.out.results.map { _run, _t, s -> s }.collect(),
        DADA2_DENOISE.out.stats.collect()
    )

    if (will_denoise && !file(params.classifier).exists()) {
        error("Classifier not found: ${params.classifier}")
    }
    ch_classifier = will_denoise ? channel.value(file(params.classifier)) : channel.empty()
    CLASSIFY(MERGE_RUNS.out.repseqs, ch_classifier)
    FILTER_TAXA(MERGE_RUNS.out.table, MERGE_RUNS.out.repseqs, CLASSIFY.out.taxonomy,
                optStr(params.exclude_taxa), optStr(params.include_taxa))
    PHYLOGENY(FILTER_TAXA.out.repseqs)

    // sample metadata: user file, or a minimal one (sample-id + run) so barplots always work
    def meta_file = params.metadata ? file(params.metadata, checkIfExists: true) : null
    if (!meta_file) {
        meta_file = file("${workflow.workDir}/sample_metadata.${workflow.sessionId}.tsv")
        meta_file.text = (['sample-id\trun'] + samples.collect { meta, _r1, _r2 -> "${meta.id}\t${meta.run}" }).join('\n') + '\n'
    }
    ch_meta = channel.value(meta_file)

    QIIME_VISUALIZE(MERGE_RUNS.out.table, FILTER_TAXA.out.table, FILTER_TAXA.out.repseqs,
                    CLASSIFY.out.taxonomy, ch_meta)
    EXPORT_TABLES(FILTER_TAXA.out.table, FILTER_TAXA.out.repseqs, CLASSIFY.out.taxonomy, PHYLOGENY.out.rooted)

    ch_div = channel.empty()
    if (params.sampling_depth > 0) {
        DIVERSITY(FILTER_TAXA.out.table, PHYLOGENY.out.rooted, ch_meta, params.sampling_depth)
        ch_div = DIVERSITY.out.core.mix(DIVERSITY.out.rarefaction, DIVERSITY.out.skipped)
        DIVERSITY.out.skipped.subscribe { f -> log.warn(f.text.trim()) }
    }

    workflow.onComplete = {
        if (workflow.success && !trunc_given && !params.auto_trunc) {
            log.info("""
                ------------------------------------------------------------------
                Stopped after QC (no truncation lengths given). Review:
                  ${workflow.outputDir}/trunc/trunc_suggestion.txt
                  ${workflow.outputDir}/qiime2/demux/*.demux.qzv   (view.qiime2.org)
                then re-run the SAME command with:  -resume --trunc_f <F> --trunc_r <R>
                ------------------------------------------------------------------""".stripIndent())
        }
    }

    publish:
    qc        = MULTIQC.out.report.mix(MULTIQC.out.data)
    trunc     = SUGGEST_TRUNC.out.report.mix(SUGGEST_TRUNC.out.tsv)
    demux     = QIIME_IMPORT.out.qza.map { _run, q -> q }.mix(QIIME_IMPORT.out.qzv)
    dada2     = DADA2_DENOISE.out.qzv.mix(MERGE_RUNS.out.stats)
    artifacts = MERGE_RUNS.out.table.mix(MERGE_RUNS.out.repseqs, CLASSIFY.out.taxonomy,
                    FILTER_TAXA.out.table, FILTER_TAXA.out.repseqs, PHYLOGENY.out.rooted)
    visuals   = QIIME_VISUALIZE.out.qzv.flatten()
    exported  = EXPORT_TABLES.out.files.flatten()
    diversity = ch_div

}

output {
    qc {
        path 'qc'
    }
    trunc {
        path 'trunc'
    }
    demux {
        path 'qiime2/demux'
    }
    dada2 {
        path 'qiime2/dada2'
    }
    artifacts {
        path 'qiime2/artifacts'
    }
    visuals {
        path 'qiime2/visualizations'
    }
    exported {
        path '.'        // files already live under exported/
    }
    diversity {
        path 'qiime2/diversity'
    }
}

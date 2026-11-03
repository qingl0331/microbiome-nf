#!/usr/bin/env nextflow
/*
 * CBI shotgun metagenome (+ optional matched metatranscriptome) pipeline
 * Replaces metagenomeHumann3/mg.master.singu.dsl2.hn3.nf and the DSL1
 * microbiomeSingularityPipeline/{metagnom,metatrans}.master.nf.
 *
 *   [merge lanes] -> [clumpify optical dedup] -> fastp -> bowtie2 host removal [-> PhiX removal]
 *   RNA only: -> SortMeRNA (index built once, cached)
 *   -> [fastq_screen QC] -> MetaPhlAn 4 (DNA; RNA without a matched DNA sample)
 *   -> HUMAnN (RNA uses the matched DNA MetaPhlAn profile) -> merged/regrouped/normalised tables
 *   -> read tracking + MultiQC
 *
 *   nextflow run shotgun.nf -profile chpc --input samples.csv
 *   nextflow run shotgun.nf -profile chpc --input samples.csv --skip_humann          # taxonomy only
 *   nextflow run shotgun.nf -profile chpc --input samples.csv --host_index /path/mm10.standard
 */

include { checkParams; optStr; readSamplesheet; readGlob; FASTQ_SCREEN; MULTIQC } from './modules/common.nf'
include {
    MERGE_LANES; CLUMPIFY; FASTP; HOST_REMOVE; SORTMERNA_INDEX; SORTMERNA; METAPHLAN; HUMANN;
    MERGE_METAPHLAN; MERGE_HUMANN; READ_TRACKING
} from './modules/shotgun.nf'

// ---- input (one of) ----
params.input            = null    // CSV: sample,fastq_1,fastq_2[,type,pair_id]; type = DNA (default) | RNA
                                  // pair_id (RNA rows) = sample name of the matched DNA library
params.reads            = null    // DNA-only shortcut glob; sample = text before first '_'

// ---- preprocessing ----
params.optical_dedup    = true
params.optical_dupedist = 12000   // 12000 NovaSeq; 2500 HiSeq 3000/4000; 40 NextSeq
params.host_index       = null    // bowtie2 index PREFIX (human default in conf/chpc.config)
params.phix_index       = null    // bowtie2 index PREFIX for PhiX spike-in removal (default in conf/chpc.config; unset = skip)
params.sortmerna_ref    = null    // rRNA fasta (smr_v4.3_default_db.fasta); required when RNA samples exist
params.db_cache         = 'db_cache' // where one-time indexes (SortMeRNA) are stored and reused

// ---- profiling ----
params.metaphlan_db     = null
params.metaphlan_index  = 'mpa_vJun23_CHOCOPhlAnSGB_202403'   // must match the HUMAnN version
params.skip_humann      = false
params.humann_nucleotide_db = null
params.humann_protein_db    = null
params.humann_utility_db    = null  // optional: enables KO / EC regrouping

// ---- QC ----
params.fastq_screen_conf = null

workflow {
    main:
    // ------------------------------------------------------------------ inputs + validation
    checkParams(['input', 'reads', 'optical_dedup', 'optical_dupedist', 'host_index', 'phix_index',
                 'sortmerna_ref', 'db_cache', 'metaphlan_db', 'metaphlan_index', 'skip_humann',
                 'humann_nucleotide_db', 'humann_protein_db', 'humann_utility_db', 'fastq_screen_conf',
                 'db_root'])   // db_root: base of the DB paths in conf/chpc.config
    if (!params.input && !params.reads) {
        error("Provide --input samplesheet.csv or --reads '<glob>'")
    }
    def samples = params.input
        ? readSamplesheet(file(params.input, checkIfExists: true), ['type', 'pair_id'])
        : readGlob(params.reads as String, [type: 'DNA', pair_id: ''])
    samples = samples.collect { meta, r1, r2 ->
        def t = (meta.type ?: 'DNA').toUpperCase()
        if (!(t in ['DNA', 'RNA'])) {
            error("Sample ${meta.id}: type must be DNA or RNA, got '${meta.type}'")
        }
        [meta + [type: t, pair_id: meta.pair_id ?: ''], r1, r2]
    }
    def dna_ids = samples.findAll { s -> s[0].type == 'DNA' }.collect { s -> s[0].id } as Set
    def n_rna = samples.count { s -> s[0].type == 'RNA' }
    samples.findAll { s -> s[0].type == 'RNA' && s[0].pair_id }.each { s ->
        if (!(s[0].pair_id in dna_ids)) {
            error("RNA sample ${s[0].id}: pair_id '${s[0].pair_id}' is not a DNA sample in the samplesheet")
        }
    }
    if (!params.host_index) {
        error("--host_index (bowtie2 prefix) is required")
    }
    def host_files = files("${params.host_index}*.bt2*")
    if (!host_files) {
        error("No bowtie2 index files found for prefix ${params.host_index}")
    }
    def phix_index = optStr(params.phix_index)   // '' = do not remove PhiX
    def phix_files = []
    if (phix_index) {
        phix_files = files("${phix_index}*.bt2*")
        if (!phix_files) {
            error("No bowtie2 index files found for --phix_index ${phix_index}")
        }
    }
    if (n_rna > 0 && !params.sortmerna_ref) {
        error("RNA samples present: --sortmerna_ref is required")
    }
    if (!params.metaphlan_db) {
        error("--metaphlan_db is required")
    }
    if (!params.skip_humann && !(params.humann_nucleotide_db && params.humann_protein_db)) {
        error("--humann_nucleotide_db and --humann_protein_db are required (or use --skip_humann)")
    }

    log.info("""
        CBI shotgun pipeline
        samples   : ${samples.size()} (${samples.size() - n_rna} DNA, ${n_rna} RNA)
        host index: ${params.host_index}
        PhiX index: ${phix_index ?: 'none (PhiX not removed)'}
        MetaPhlAn : ${params.metaphlan_index}
        HUMAnN    : ${params.skip_humann ? 'skipped' : 'yes'}${params.humann_utility_db ? ' (+ KO/EC regroup)' : ''}
        """.stripIndent())

    // sample -> type table for read tracking
    def type_file = file("${workflow.workDir}/sample_types.${workflow.sessionId}.tsv")
    type_file.text = (['sample\ttype'] + samples.collect { s -> "${s[0].id}\t${s[0].type}" }).join('\n') + '\n'

    // ------------------------------------------------------------------ preprocessing
    ch_in = channel.fromList(samples).branch { _meta, r1s, _r2s ->
        multi: r1s.size() > 1
        single: true
    }
    MERGE_LANES(ch_in.multi)
    ch_raw = ch_in.single
        .map { meta, r1s, r2s -> [meta, r1s[0], r2s[0]] }
        .mix(MERGE_LANES.out.reads)

    ch_dedup = ch_raw
    ch_clump_log = channel.empty()
    if (params.optical_dedup) {
        CLUMPIFY(ch_raw, params.optical_dupedist)
        ch_dedup = CLUMPIFY.out.reads
        ch_clump_log = CLUMPIFY.out.log
    }

    FASTP(ch_dedup)
    HOST_REMOVE(
        FASTP.out.reads,
        host_files, file(params.host_index).name,
        phix_files, phix_index ? file(phix_index).name : ''
    )

    ch_clean_split = HOST_REMOVE.out.reads.branch { meta, _r1, _r2 ->
        rna: meta.type == 'RNA'
        dna: true
    }

    ch_rna_clean = channel.empty()
    ch_rrna_counts = channel.empty()
    ch_rrna_logs = channel.empty()
    if (n_rna > 0) {
        def smr_ref = file(params.sortmerna_ref, checkIfExists: true)
        SORTMERNA_INDEX(smr_ref, params.db_cache)
        SORTMERNA(ch_clean_split.rna, smr_ref, SORTMERNA_INDEX.out.idx)
        ch_rna_clean = SORTMERNA.out.reads
        ch_rrna_counts = SORTMERNA.out.counts
        ch_rrna_logs = SORTMERNA.out.log
    }
    ch_final = ch_clean_split.dna.mix(ch_rna_clean)

    ch_screen = channel.empty()
    if (params.fastq_screen_conf) {
        FASTQ_SCREEN(ch_final.map { meta, r1, _r2 -> [meta, r1] }, file(params.fastq_screen_conf, checkIfExists: true))
        ch_screen = FASTQ_SCREEN.out.txt
    }

    // ------------------------------------------------------------------ taxonomy
    // MetaPhlAn on every DNA sample, and on RNA samples that have no matched DNA library
    ch_for_mpa = ch_final.filter { meta, _r1, _r2 -> meta.type == 'DNA' || !meta.pair_id }
    METAPHLAN(ch_for_mpa, file(params.metaphlan_db, checkIfExists: true), params.metaphlan_index)

    MERGE_METAPHLAN(
        METAPHLAN.out.profile
            .map { meta, p -> [meta.type, p] }
            .groupTuple()
    )

    // ------------------------------------------------------------------ function
    ch_humann_tables = channel.empty()
    ch_humann_logs = channel.empty()
    if (!params.skip_humann) {
        ch_profiles = METAPHLAN.out.profile.map { meta, p -> [meta.id, p] }
        // key each sample by the sample whose profile it should use
        ch_humann_in = ch_final
            .map { meta, r1, r2 -> [(meta.type == 'RNA' && meta.pair_id) ? meta.pair_id : meta.id, meta, r1, r2] }
            .combine(ch_profiles, by: 0)
            .map { _key, meta, r1, r2, p -> [meta, r1, r2, p] }

        HUMANN(ch_humann_in, file(params.humann_nucleotide_db, checkIfExists: true), file(params.humann_protein_db, checkIfExists: true))

        def ch_util = params.humann_utility_db
            ? channel.value(file(params.humann_utility_db, checkIfExists: true))
            : channel.value(file("${projectDir}/assets/NO_UTILITY_DB"))
        MERGE_HUMANN(
            HUMANN.out.genefamilies.map { meta, f -> [meta.type, f] }.groupTuple()
                .join(HUMANN.out.pathabundance.map { meta, f -> [meta.type, f] }.groupTuple())
                .join(HUMANN.out.pathcoverage.map { meta, f -> [meta.type, f] }.groupTuple(), remainder: true)
                .map { type, gf, pa, pc -> [type, gf, pa, pc ?: []] },
            ch_util
        )
        ch_humann_tables = MERGE_HUMANN.out.tables.flatten()
        ch_humann_logs = HUMANN.out.log
    }

    // ------------------------------------------------------------------ reporting
    READ_TRACKING(
        FASTP.out.json.collect(),
        HOST_REMOVE.out.counts.mix(ch_rrna_counts).collect(),
        type_file
    )
    MULTIQC(
        FASTP.out.json
            .mix(HOST_REMOVE.out.log, ch_rrna_logs, ch_screen, ch_clump_log)
            .collect(),
        'CBI shotgun QC'
    )

    publish:
    qc            = MULTIQC.out.report.mix(MULTIQC.out.data, READ_TRACKING.out.tsv, FASTP.out.html)
    metaphlan     = MERGE_METAPHLAN.out.tables.flatten()
    metaphlan_raw = METAPHLAN.out.profile.map { _meta, p -> p }.mix(METAPHLAN.out.bt2out)
    humann        = ch_humann_tables
    humann_logs   = ch_humann_logs
}

output {
    qc {
        path 'qc'
    }
    metaphlan {
        path 'metaphlan'
    }
    metaphlan_raw {
        path 'metaphlan/per_sample'
    }
    humann {
        path '.'        // humann_DNA/ and humann_RNA/
    }
    humann_logs {
        path 'humann_logs'
    }
}

/*
 * Shotgun metagenome / metatranscriptome processes
 * (replaces metagenomeHumann3/ and the DSL1 microbiomeSingularityPipeline/)
 */

/*
 * Concatenate lanes — only runs for samples that actually have >1 lane.
 */
process MERGE_LANES {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(r1s, stageAs: 'in_R1/*'), path(r2s, stageAs: 'in_R2/*')

    output:
    tuple val(meta), path("${meta.id}_R1.merged.fastq.gz"), path("${meta.id}_R2.merged.fastq.gz"), emit: reads

    script:
    """
    cat in_R1/* > ${meta.id}_R1.merged.fastq.gz
    cat in_R2/* > ${meta.id}_R2.merged.fastq.gz
    """

    stub:
    """
    touch ${meta.id}_R1.merged.fastq.gz ${meta.id}_R2.merged.fastq.gz
    """
}

/*
 * Optical-duplicate removal (patterned flow cells), as before: dupedist=12000 for NovaSeq.
 */
process CLUMPIFY {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(r1), path(r2)
    val dupedist

    output:
    tuple val(meta), path("${meta.id}_R1.clump.fastq.gz"), path("${meta.id}_R2.clump.fastq.gz"), emit: reads
    path "${meta.id}.clumpify.log", emit: log

    script:
    def mem = Math.max(2, (task.memory.toGiga() * 0.85) as int)
    """
    clumpify.sh -Xmx${mem}g threads=${task.cpus} \\
        in1=${r1} in2=${r2} \\
        out1=${meta.id}_R1.clump.fastq.gz out2=${meta.id}_R2.clump.fastq.gz \\
        dedupe=t optical=t dupedist=${dupedist} 2> ${meta.id}.clumpify.log
    """

    stub:
    """
    touch ${meta.id}_R1.clump.fastq.gz ${meta.id}_R2.clump.fastq.gz ${meta.id}.clumpify.log
    """
}

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.fastp.fastq.gz"), path("${meta.id}_R2.fastp.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    def args = task.ext.args ?: ''
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${meta.id}_R1.fastp.fastq.gz -O ${meta.id}_R2.fastp.fastq.gz \\
        --thread ${Math.min(task.cpus as int, 16)} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_R1.fastp.fastq.gz ${meta.id}_R2.fastp.fastq.gz ${meta.id}.fastp.html
    echo '{"summary":{"before_filtering":{"total_reads":2000},"after_filtering":{"total_reads":1800}}}' > ${meta.id}.fastp.json
    """
}

/*
 * Host read removal with bowtie2 against any host index (human, mouse, ...), then (optional)
 * a second, cheap bowtie2 pass against PhiX — the old kneaddata step removed hg38 + PhiX.
 * In each pass a pair is kept only if BOTH mates are unmapped (samtools fastq -f 12 -F 256).
 */
process HOST_REMOVE {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)
    path index_files, stageAs: 'host_index/*'
    val index_name
    path phix_files, stageAs: 'phix_index/*'   // [] when --phix_index is not set
    val phix_name                               // '' = skip the PhiX pass

    output:
    tuple val(meta), path("${meta.id}_R1.host_removed.fastq.gz"), path("${meta.id}_R2.host_removed.fastq.gz"), emit: reads
    path "*.bowtie2.log", emit: log             // <id>.bowtie2.log (host) [+ <id>_phix.bowtie2.log]
    path "${meta.id}.host_removed.counts.tsv", emit: counts

    script:
    def args = task.ext.args ?: ''
    def host_out = phix_name ? 'host_clean' : "${meta.id}"
    def phix_pass = phix_name ? """
    bowtie2 -p ${task.cpus} ${args} -x phix_index/${phix_name} -1 host_clean_R1.host_removed.fastq.gz -2 host_clean_R2.host_removed.fastq.gz \\
        2> ${meta.id}_phix.bowtie2.log \\
      | samtools fastq -@ 2 -f 12 -F 256 \\
            -1 ${meta.id}_R1.host_removed.fastq.gz -2 ${meta.id}_R2.host_removed.fastq.gz \\
            -0 /dev/null -s /dev/null -
    rm host_clean_R1.host_removed.fastq.gz host_clean_R2.host_removed.fastq.gz
    n=\$(( \$(pigz -dc ${meta.id}_R1.host_removed.fastq.gz | wc -l) / 4 ))
    printf '%s\\tphix_removed\\t%s\\n' "${meta.id}" "\$n" >> ${meta.id}.host_removed.counts.tsv
    """ : ''
    """
    set -o pipefail   # a bowtie2 failure must not silently yield empty "clean" reads
    bowtie2 -p ${task.cpus} ${args} -x host_index/${index_name} -1 ${r1} -2 ${r2} 2> ${meta.id}.bowtie2.log \\
      | samtools fastq -@ 2 -f 12 -F 256 \\
            -1 ${host_out}_R1.host_removed.fastq.gz -2 ${host_out}_R2.host_removed.fastq.gz \\
            -0 /dev/null -s /dev/null -
    n=\$(( \$(pigz -dc ${host_out}_R1.host_removed.fastq.gz | wc -l) / 4 ))
    printf '%s\\thost_removed\\t%s\\n' "${meta.id}" "\$n" > ${meta.id}.host_removed.counts.tsv
    ${phix_pass}
    """

    stub:
    def phix_stub = phix_name ? "touch ${meta.id}_phix.bowtie2.log; printf '%s\\tphix_removed\\t1490\\n' \"${meta.id}\" >> ${meta.id}.host_removed.counts.tsv" : ''
    """
    touch ${meta.id}_R1.host_removed.fastq.gz ${meta.id}_R2.host_removed.fastq.gz ${meta.id}.bowtie2.log
    printf '%s\\thost_removed\\t1500\\n' "${meta.id}" > ${meta.id}.host_removed.counts.tsv
    ${phix_stub}
    """
}

/*
 * Build the SortMeRNA index once and keep it (storeDir) — the old pipeline rebuilt it per sample.
 */
process SORTMERNA_INDEX {
    label 'process_medium'
    storeDir "${cache_dir}/sortmerna_${ref.baseName}"

    input:
    path ref
    val cache_dir

    output:
    path "idx", emit: idx

    script:
    """
    sortmerna --ref ${ref} --index 1 --workdir . --idx-dir idx --threads ${task.cpus}
    """

    stub:
    """
    mkdir idx
    """
}

/*
 * rRNA depletion for metatranscriptomes. --paired_in: a pair counts as rRNA if either mate hits,
 * so the "other" output keeps only clean pairs (no PEAR stitching needed any more).
 */
process SORTMERNA {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)
    path ref
    path idx

    output:
    tuple val(meta), path("${meta.id}_R1.non_rRNA.fastq.gz"), path("${meta.id}_R2.non_rRNA.fastq.gz"), emit: reads
    path "${meta.id}.sortmerna.log", emit: log
    path "${meta.id}.non_rRNA.counts.tsv", emit: counts

    script:
    """
    sortmerna --ref ${ref} --index 0 --idx-dir ${idx} \\
        --reads ${r1} --reads ${r2} \\
        --workdir . --threads ${task.cpus} \\
        --paired_in --out2 --fastx \\
        --aligned rRNA_reads --other non_rRNA_reads

    mv non_rRNA_reads_fwd.f*q.gz ${meta.id}_R1.non_rRNA.fastq.gz
    mv non_rRNA_reads_rev.f*q.gz ${meta.id}_R2.non_rRNA.fastq.gz
    mv rRNA_reads.log ${meta.id}.sortmerna.log
    rm -rf kvdb readb rRNA_reads_*.f*q.gz

    n=\$(( \$(gzip -dc ${meta.id}_R1.non_rRNA.fastq.gz | wc -l) / 4 ))
    printf '%s\\tnon_rRNA\\t%s\\n' "${meta.id}" "\$n" > ${meta.id}.non_rRNA.counts.tsv
    """

    stub:
    """
    touch ${meta.id}_R1.non_rRNA.fastq.gz ${meta.id}_R2.non_rRNA.fastq.gz ${meta.id}.sortmerna.log
    printf '%s\\tnon_rRNA\\t900\\n' "${meta.id}" > ${meta.id}.non_rRNA.counts.tsv
    """
}

/*
 * MetaPhlAn run on its own (not inside HUMAnN): both mates are used, the profile is reusable
 * for taxonomy-only projects, and the bowtie2out allows cheap re-profiling with other options.
 */
process METAPHLAN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)
    path db
    val index

    output:
    tuple val(meta), path("${meta.id}_metaphlan_profile.tsv"), emit: profile
    path "${meta.id}.bowtie2out.bz2", emit: bt2out

    script:
    def args = task.ext.args ?: ''
    """
    metaphlan ${r1},${r2} \\
        --input_type fastq \\
        --nproc ${task.cpus} \\
        --bowtie2db ${db} \\
        --index ${index} \\
        --bowtie2out ${meta.id}.bowtie2out.bz2 \\
        --sample_id_key sample_id --sample_id ${meta.id} \\
        ${args} \\
        -o ${meta.id}_metaphlan_profile.tsv
    """

    stub:
    """
    printf '#mpa_vJun23_CHOCOPhlAnSGB_202403\\n#clade_name\\tNCBI_tax_id\\trelative_abundance\\tadditional_species\\n' > ${meta.id}_metaphlan_profile.tsv
    touch ${meta.id}.bowtie2out.bz2
    """
}

/*
 * HUMAnN with the pre-computed MetaPhlAn profile (DNA: own profile; RNA: matched DNA profile).
 * Output file globs cover both HUMAnN 3.9 and 4.x naming.
 */
process HUMANN {
    tag "${meta.id}"
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(r1), path(r2), path(profile)
    path nucleotide_db
    path protein_db

    output:
    tuple val(meta), path("${meta.id}*genefamilies.tsv.gz"), emit: genefamilies
    tuple val(meta), path("${meta.id}*pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("${meta.id}*pathcoverage.tsv.gz"), emit: pathcoverage, optional: true
    path "${meta.id}.humann.log", emit: log

    script:
    def args = task.ext.args ?: ''
    """
    cat ${r1} ${r2} > ${meta.id}.fastq.gz

    humann \\
        --input ${meta.id}.fastq.gz \\
        --output humann_out \\
        --output-basename ${meta.id} \\
        --taxonomic-profile ${profile} \\
        --nucleotide-database ${nucleotide_db} \\
        --protein-database ${protein_db} \\
        --threads ${task.cpus} \\
        --o-log ${meta.id}.humann.log \\
        ${args}

    for f in humann_out/*.tsv; do gzip -c "\$f" > "\$(basename \$f).gz"; done
    rm -rf humann_out ${meta.id}.fastq.gz
    """

    stub:
    """
    for t in genefamilies pathabundance pathcoverage; do echo "" | gzip > ${meta.id}_\${t}.tsv.gz; done
    touch ${meta.id}.humann.log
    """
}

process MERGE_METAPHLAN {
    tag "${type}"
    label 'process_single'

    input:
    tuple val(type), path(profiles)

    output:
    path "metaphlan_${type}_*.tsv", emit: tables

    script:
    """
    merge_metaphlan_tables.py ${profiles} > metaphlan_${type}_merged.tsv
    # species-level only (SGB/t__ rows removed), handy for R
    grep -E '^(clade_name|#)|s__' metaphlan_${type}_merged.tsv | grep -v 't__' > metaphlan_${type}_species.tsv
    """

    stub:
    """
    touch metaphlan_${type}_merged.tsv metaphlan_${type}_species.tsv
    """
}

/*
 * Join per-sample HUMAnN tables, regroup UniRef90 -> KEGG KO / EC (if utility DB given),
 * normalise (gene families/KO/EC -> CPM, pathways -> relab + CPM), split strat/unstrat.
 */
process MERGE_HUMANN {
    tag "${type}"
    label 'process_medium'

    input:
    tuple val(type), path(gf, stageAs: 'gf/*'), path(pa, stageAs: 'pa/*'), path(pc, stageAs: 'pc/*')
    path utility_db

    output:
    path "humann_${type}/*.tsv", emit: tables

    script:
    def util = utility_db.name != 'NO_UTILITY_DB' ? utility_db : ''
    """
    mkdir -p j_gf j_pa j_pc humann_${type}
    for f in gf/*.gz; do gzip -dc "\$f" > j_gf/\$(basename "\${f%.gz}"); done
    for f in pa/*.gz; do gzip -dc "\$f" > j_pa/\$(basename "\${f%.gz}"); done
    for f in pc/*.gz; do [ -e "\$f" ] && gzip -dc "\$f" > j_pc/\$(basename "\${f%.gz}"); done || true

    out=humann_${type}
    humann_join_tables -i j_gf -o \$out/genefamilies_rpk.tsv --file_name genefamilies
    humann_join_tables -i j_pa -o \$out/pathabundance.tsv --file_name pathabundance
    if ls j_pc/*.tsv >/dev/null 2>&1; then
        humann_join_tables -i j_pc -o \$out/pathcoverage.tsv --file_name pathcoverage
    fi

    humann_renorm_table -i \$out/genefamilies_rpk.tsv -o \$out/genefamilies_cpm.tsv --units cpm --update-snames
    humann_renorm_table -i \$out/pathabundance.tsv -o \$out/pathabundance_relab.tsv --units relab --update-snames
    humann_renorm_table -i \$out/pathabundance.tsv -o \$out/pathabundance_cpm.tsv --units cpm --update-snames

    util="${util}"
    if [ -n "\$util" ]; then
        for g in ko level4ec; do
            if [ -e "\$util/map_\${g}_uniref90.txt.gz" ]; then
                humann_regroup_table -i \$out/genefamilies_rpk.tsv -c "\$util/map_\${g}_uniref90.txt.gz" -o \$g.rpk.tsv
                humann_renorm_table -i \$g.rpk.tsv -o \$g.cpm.tsv --units cpm --update-snames
                n=\$g; [ "\$g" = level4ec ] && n=ec      # HUMAnN 3 ships map_ec_name, not map_level4ec_name
                if [ -e "\$util/map_\${n}_name.txt.gz" ]; then
                    humann_rename_table -i \$g.cpm.tsv -c "\$util/map_\${n}_name.txt.gz" -o \$out/\${g}_cpm.tsv
                else
                    mv \$g.cpm.tsv \$out/\${g}_cpm.tsv
                fi
            fi
        done
    fi

    mkdir -p split_tmp
    for t in \$out/*.tsv; do
        humann_split_stratified_table -i "\$t" -o split_tmp
    done
    mv split_tmp/*.tsv \$out/
    """

    stub:
    """
    mkdir -p humann_${type}
    touch humann_${type}/genefamilies_cpm.tsv humann_${type}/pathabundance_relab.tsv
    """
}

/*
 * Per-sample read accounting across every step (fastp in/out, host removal, rRNA removal).
 */
process READ_TRACKING {
    label 'process_single'

    input:
    path jsons, stageAs: 'fastp/*'
    path counts, stageAs: 'counts/*'
    path samplesheet_tsv

    output:
    path "read_tracking.tsv", emit: tsv

    script:
    """
    read_tracking.py --fastp fastp/*.json --counts counts/* --samples ${samplesheet_tsv} > read_tracking.tsv
    """

    stub:
    """
    touch read_tracking.tsv
    """
}

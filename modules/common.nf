/*
 * Shared helpers + processes used by both amplicon.nf and shotgun.nf
 */

/*
 * Classic params accept any name, so a typo (e.g. --trunclen) would be silently ignored.
 * Fail instead on any --option on the command line that the entry script does not define.
 */
def checkParams(List<String> known) {
    def given = (workflow.commandLine =~ /(?:^|\s)--([A-Za-z][\w-]*)/).collect { m -> m[1] }
    def unknown = given.findAll { p -> !(p in known) }
    if (unknown) {
        error("Unknown parameter(s): ${unknown.collect { p -> "--${p}" }.join(', ')} -- check the spelling against the params list in ${workflow.scriptName}")
    }
}

/*
 * An empty value on the command line (--phix_index '') reaches classic params as Boolean true.
 * Use for string params where '' means "off".
 */
def optStr(value) {
    return (value instanceof Boolean || value == null) ? '' : value.toString()
}

/*
 * Read a samplesheet eagerly (so it can be validated before anything is submitted).
 * Required columns: sample, fastq_1, fastq_2. Optional columns are passed through.
 * Several rows with the same `sample` = several lanes/runs of one library -> concatenated.
 * Returns: List of [meta(Map), List<Path> r1s, List<Path> r2s]
 */
def readSamplesheet(Path sheet, List<String> optionalCols) {
    def rows = sheet.splitCsv(header: true, strip: true)
    if (!rows) {
        error("Samplesheet ${sheet} is empty")
    }
    def required = ['sample', 'fastq_1', 'fastq_2']
    def missing = required.findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet ${sheet} is missing column(s): ${missing.join(', ')}")
    }
    def bySample = [:]
    rows.each { row ->
        def id = row.sample
        if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
            error("Sample name '${id}' must start with a letter/number and contain only letters, numbers, '.', '_' or '-'")
        }
        def r1 = file(row.fastq_1, checkIfExists: true)
        def r2 = file(row.fastq_2, checkIfExists: true)
        def extra = optionalCols.collectEntries { col -> [(col): (row[col] ?: '')] }
        if (!bySample.containsKey(id)) {
            bySample[id] = [meta: [id: id] + extra, r1: [], r2: []]
        }
        else if (bySample[id].meta.subMap(optionalCols) != extra) {
            error("Sample '${id}' appears on several rows with different ${optionalCols.join('/')} values")
        }
        bySample[id].r1 << r1
        bySample[id].r2 << r2
    }
    return bySample.values().collect { rec -> [rec.meta, rec.r1, rec.r2] }
}

/*
 * Glob mode (quick runs, no samplesheet): sample ID = text before the first '_' of the
 * file name (GNomEx style, e.g. 15215X1_..._R1_001.fastq.gz -> 15215X1). Files sharing a
 * prefix (lanes) are merged, matching the old manifest behaviour.
 */
def readGlob(String pattern, Map extraMeta) {
    def pairs = [:]
    files(pattern, checkIfExists: true).each { f ->
        def name = f.name
        def m = name =~ /^(.+?)_.*?_?R([12])(_\d{3})?\.f(ast)?q\.gz$/
        if (!m.find()) {
            error("Cannot tell R1/R2 apart for ${name}; use a samplesheet (--input) instead")
        }
        def id = name.tokenize('_')[0]
        def mate = m.group(2)
        if (!pairs.containsKey(id)) {
            pairs[id] = [r1: [], r2: []]
        }
        pairs[id]["r${mate}"] << f
    }
    if (!pairs) {
        error("No FASTQ files match ${pattern}")
    }
    return pairs.collect { id, rec ->
        if (rec.r1.size() != rec.r2.size()) {
            error("Sample ${id}: ${rec.r1.size()} R1 files but ${rec.r2.size()} R2 files")
        }
        [[id: id] + extraMeta, rec.r1.sort(), rec.r2.sort()]
    }
}

process FASTQ_SCREEN {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path conf

    output:
    path "*_screen.txt", emit: txt
    path "*_screen.html", emit: html

    script:
    def args = task.ext.args ?: ''
    """
    fastq_screen --conf ${conf} --threads ${task.cpus} --aligner bowtie2 ${args} ${reads}
    """

    stub:
    """
    touch ${meta.id}_screen.txt ${meta.id}_screen.html
    """
}

process MULTIQC {
    label 'process_single'

    input:
    path files, stageAs: 'mqc_in/??/*'
    val title

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc --force --title "${title}" --filename multiqc_report.html mqc_in
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data
    """
}

#!/usr/bin/env nextflow
/*
 * Train a naive Bayes classifier for the QIIME 2 version pinned in conf/containers.config.
 * Needed because classifiers are tied to the scikit-learn version, and since QIIME 2 2026.4
 * pre-trained classifiers are only distributed externally (SILVA site, GTDB on Zenodo) —
 * there is no Greengenes2 classifier for 2026.x, so GG2 must be trained here.
 *
 * Greengenes2 example (files from the GG2 2024.09 release):
 *   nextflow run train_classifier.nf -profile chpc \
 *     --ref_seqs 2024.09.backbone.full-length.fna.qza --ref_tax 2024.09.backbone.tax.qza \
 *     --name gg2_2024.09_full_length -output-dir $DB/qiime2_2026.7/classifiers
 *
 * Full-length training is the current QIIME 2 recommendation; add --primer_f/--primer_r
 * to extract the amplicon region instead.
 */

include { checkParams; optStr } from './modules/common.nf'

params.ref_seqs = null    // .qza FeatureData[Sequence] or .fasta            (required)
params.ref_tax  = null    // .qza FeatureData[Taxonomy] or headerless 2-column .tsv/.txt (required)
params.name     = null    // classifier name -> <name>.nb.qza               (required)
params.primer_f = ''
params.primer_r = ''

process IMPORT_REFERENCE {
    label 'process_low'

    input:
    path seqs
    path tax

    output:
    tuple path("ref_seqs.qza"), path("ref_tax.qza"), emit: ref

    script:
    """
    export TMPDIR="\$PWD/tmp" HOME="\$PWD"; mkdir -p "\$TMPDIR"
    case "${seqs}" in
      *.qza) cp -L ${seqs} ref_seqs.qza ;;
      *) qiime tools import --type 'FeatureData[Sequence]' --input-path ${seqs} --output-path ref_seqs.qza ;;
    esac
    case "${tax}" in
      *.qza) cp -L ${tax} ref_tax.qza ;;
      *) qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat \\
             --input-path ${tax} --output-path ref_tax.qza ;;
    esac
    """

    stub:
    """
    touch ref_seqs.qza ref_tax.qza
    """
}

process TRAIN_CLASSIFIER {
    label 'process_single'
    label 'process_high_memory'
    label 'process_long'

    input:
    tuple path(seqs), path(tax)
    val name
    val primer_f
    val primer_r

    output:
    path "${name}.nb.qza", emit: classifier
    path "${name}.versions.txt", emit: versions

    script:
    def extract = primer_f && primer_r
    """
    export TMPDIR="\$PWD/tmp" HOME="\$PWD"; mkdir -p "\$TMPDIR"
    reads=${seqs}
    if ${extract}; then
        qiime feature-classifier extract-reads --i-sequences ${seqs} \\
            --p-f-primer ${primer_f} --p-r-primer ${primer_r} --p-n-jobs ${task.cpus} \\
            --o-reads region_seqs.qza
        reads=region_seqs.qza
    fi
    qiime feature-classifier fit-classifier-naive-bayes \\
        --i-reference-reads \$reads \\
        --i-reference-taxonomy ${tax} \\
        --o-classifier ${name}.nb.qza
    { qiime info; python -c 'import sklearn; print("scikit-learn", sklearn.__version__)'; } > ${name}.versions.txt
    """

    stub:
    """
    touch ${name}.nb.qza ${name}.versions.txt
    """
}

workflow {
    main:
    checkParams(['ref_seqs', 'ref_tax', 'name', 'primer_f', 'primer_r'])
    if (!(params.ref_seqs && params.ref_tax && params.name)) {
        error("--ref_seqs, --ref_tax and --name are required")
    }
    IMPORT_REFERENCE(file(params.ref_seqs, checkIfExists: true), file(params.ref_tax, checkIfExists: true))
    TRAIN_CLASSIFIER(IMPORT_REFERENCE.out.ref, params.name, optStr(params.primer_f), optStr(params.primer_r))

    publish:
    classifier = TRAIN_CLASSIFIER.out.classifier.mix(TRAIN_CLASSIFIER.out.versions)
}

output {
    classifier {
        path '.'
    }
}

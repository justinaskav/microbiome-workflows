// Length filtering removed - not needed for UNITE reference sequences
// ITS regions are naturally variable in length and UNITE sequences are already curated

process DEREPLICATE_UNITE {
    tag "$meta.id"
    label 'process_high'
    cpus 2

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val derep_mode_str

    output:
    tuple val(meta), path("${meta.id}-seqs-derep.qza"), path("${meta.id}-tax-derep.qza"), emit: dereplicated_data

    script:
    """
    qiime rescript dereplicate \
        --i-sequences "${sequences_qza}" \
        --i-taxa "${taxonomy_qza}" \
        --p-mode '${derep_mode_str}' \
        --p-threads ${task.cpus} \
        --o-dereplicated-sequences "${meta.id}-seqs-derep.qza" \
        --o-dereplicated-taxa "${meta.id}-tax-derep.qza" \
        --verbose
    """
}

process EXTRACT_ITS_REGION {
    tag "$meta.id (ITS Region)"
    label 'process_medium'

    conda "${params.qiime2_conda_env}"
    cpus 4

    input:
    tuple val(meta), path(sequences_qza)
    val its_region

    output:
    tuple val(meta), path("${meta.id}-its-extracted.qza"), emit: extracted_sequences

    script:
    def region_param = its_region.toLowerCase() == 'its1' ? 'ITS1' : 'ITS2'
    """
    qiime feature-classifier extract-reads \
        --i-sequences "${sequences_qza}" \
        --p-f-primer 'CTTGGTCATTTAGAGGAAGTAA' \
        --p-r-primer 'GCTGCGTTCTTCATCGATGC' \
        --p-n-jobs ${task.cpus} \
        --p-read-orientation 'forward' \
        --o-reads "${meta.id}-its-extracted.qza" \
        --verbose
    """
}

process CREATE_UNITE_CLASSIFIER {
    tag "$meta.id (UNITE Classifier)"
    label 'process_high'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(reference_reads_qza), path(reference_taxonomy_qza)

    output:
    path("${meta.id}-classifier.qza"), emit: classifier

    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "${reference_reads_qza}" \
        --i-reference-taxonomy "${reference_taxonomy_qza}" \
        --p-feat-ext--ngram-range '[8, 8]' \
        --o-classifier "${meta.id}-classifier.qza" \
        --verbose
    """
}

include { DOWNLOAD_UNITE } from './download_unite_wf.nf'

workflow PRECLASSIFY_UNITE {
    take:
    meta
    unite_version
    unite_taxon_group
    unite_cluster_id
    include_singletons
    its_region
    derep_mode

    main:
    // Download UNITE data
    unite_data_ch = DOWNLOAD_UNITE(
        meta,
        unite_version,
        unite_taxon_group,
        unite_cluster_id,
        include_singletons
    )

    // Extract ITS region (if specified)
    if (its_region && its_region != 'all') {
        extracted_seqs_ch = EXTRACT_ITS_REGION(
            unite_data_ch.sequences,
            its_region
        ).extracted_sequences
        
        processing_seqs_ch = extracted_seqs_ch
    } else {
        processing_seqs_ch = unite_data_ch.sequences
    }

    // Skip length filtering - UNITE sequences are already curated
    // ITS regions are naturally variable in length, filtering could remove valid sequences
    
    // Prepare input for dereplication: [meta, seqs.qza, tax.qza]
    derep_input_ch = processing_seqs_ch.join(unite_data_ch.taxonomy, by: 0)

    // Dereplicate sequences and taxonomy
    dereplicated_ch = DEREPLICATE_UNITE(
        derep_input_ch,
        derep_mode
    ).dereplicated_data

    // Create classifier
    classifier = CREATE_UNITE_CLASSIFIER(dereplicated_ch).classifier

    emit:
    classifier
}
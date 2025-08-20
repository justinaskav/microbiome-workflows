process DOWNLOAD_UNITE_DATA {
    tag "$meta.id"
    label 'process_low'

    conda "${params.qiime2_conda_env}"

    input:
    val meta
    val unite_version
    val unite_taxon_group
    val unite_cluster_id
    val include_singletons

    output:
    tuple val(meta), path("${meta.id}-seqs.qza"), path("${meta.id}-tax.qza"), emit: unite_data

    script:
    def singletons_flag = include_singletons ? '--p-singletons' : '--p-no-singletons'
    """
    qiime rescript get-unite-data \
        --p-version '${unite_version}' \
        --p-taxon-group '${unite_taxon_group}' \
        --p-cluster-id '${unite_cluster_id}' \
        ${singletons_flag} \
        --o-sequences "${meta.id}-seqs.qza" \
        --o-taxonomy "${meta.id}-tax.qza" \
        --verbose
    """
}

workflow DOWNLOAD_UNITE {
    take:
    meta
    unite_version
    unite_taxon_group
    unite_cluster_id
    include_singletons

    main:
    unite_data_ch = DOWNLOAD_UNITE_DATA(
        meta,
        unite_version,
        unite_taxon_group,
        unite_cluster_id,
        include_singletons
    ).unite_data

    // Split the combined output into separate channels
    sequences_ch = unite_data_ch.map { meta, seqs, tax -> tuple(meta, seqs) }
    taxonomy_ch = unite_data_ch.map { meta, seqs, tax -> tuple(meta, tax) }

    emit:
    sequences = sequences_ch
    taxonomy = taxonomy_ch
}
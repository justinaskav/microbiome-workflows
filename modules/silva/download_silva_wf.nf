// Download and initial processing of SILVA data
process GET_SILVA_DATA {
    tag "$meta.id"

    conda "${params.qiime2_conda_env}"

    input:
    val meta
    val silva_version
    val silva_target
    val include_species_labels

    output:
    tuple val(meta), path("silva-${silva_version}-${silva_target}-rna-seqs.qza"), path("silva-${silva_version}-${silva_target}-tax.qza"), emit: silva_data

    script:
    def species_flag = include_species_labels ? "--p-include-species-labels" : ""
    """
    qiime rescript get-silva-data \
        --p-version '${silva_version}' \
        --p-target '${silva_target}' \
        ${species_flag} \
        --o-silva-sequences silva-${silva_version}-${silva_target}-rna-seqs.qza \
        --o-silva-taxonomy silva-${silva_version}-${silva_target}-tax.qza \
        --verbose
    """
}

workflow DOWNLOAD_SILVA {
    take:
    meta
    silva_version
    silva_target
    include_species_labels

    main:
    // Get SILVA data and create a channel with RNA sequences and taxonomy
    raw_silva_ch = GET_SILVA_DATA(
        meta,
        silva_version,
        silva_target,
        include_species_labels
    ).silva_data
    
    emit:
    raw_silva_ch
}

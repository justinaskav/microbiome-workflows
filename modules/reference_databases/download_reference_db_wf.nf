// Unified reference database download workflow
// Supports both SILVA and UNITE databases with parameterized configuration

process DOWNLOAD_SILVA_DATA {
    tag "$meta.id - SILVA v${version}"
    label 'process_low'

    conda "${params.qiime2_conda_env}"

    input:
    val meta
    val version
    val target
    val include_species_labels

    output:
    tuple val(meta), path("silva-${version}-${target}-rna-seqs.qza"), path("silva-${version}-${target}-tax.qza"), emit: reference_data

    script:
    def species_flag = include_species_labels ? "--p-include-species-labels" : ""
    """
    qiime rescript get-silva-data \
        --p-version '${version}' \
        --p-target '${target}' \
        ${species_flag} \
        --o-silva-sequences silva-${version}-${target}-rna-seqs.qza \
        --o-silva-taxonomy silva-${version}-${target}-tax.qza \
        --verbose
    """
}

process DOWNLOAD_UNITE_DATA {
    tag "$meta.id - UNITE v${version}"
    label 'process_low'

    conda "${params.qiime2_conda_env}"

    input:
    val meta
    val version
    val taxon_group
    val cluster_id
    val include_singletons

    output:
    tuple val(meta), path("${meta.id}-seqs.qza"), path("${meta.id}-tax.qza"), emit: reference_data

    script:
    def singletons_flag = include_singletons ? '--p-singletons' : '--p-no-singletons'
    """
    qiime rescript get-unite-data \
        --p-version '${version}' \
        --p-taxon-group '${taxon_group}' \
        --p-cluster-id '${cluster_id}' \
        ${singletons_flag} \
        --o-sequences "${meta.id}-seqs.qza" \
        --o-taxonomy "${meta.id}-tax.qza" \
        --verbose
    """
}

process DOWNLOAD_DADA2_SILVA_REFERENCES {
    tag "DADA2 SILVA references"
    label 'process_low'
    
    output:
    path("silva_nr99_v138.2_toGenus_trainset.fa.gz"), emit: trainset
    path("silva_v138.2_assignSpecies.fa.gz"), emit: species
    
    script:
    """
    curl -L -o silva_nr99_v138.2_toGenus_trainset.fa.gz "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"
    curl -L -o silva_v138.2_assignSpecies.fa.gz "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1"
    """
}

workflow DOWNLOAD_REFERENCE_DATABASE {
    take:
    meta
    db_type              // 'silva', 'unite', or 'dada2_silva'
    db_config            // Map containing database-specific parameters

    main:
    if (db_type == 'silva') {
        reference_data_ch = DOWNLOAD_SILVA_DATA(
            meta,
            db_config.version,
            db_config.target,
            db_config.include_species_labels ?: false
        ).reference_data
        
    } else if (db_type == 'unite') {
        reference_data_ch = DOWNLOAD_UNITE_DATA(
            meta,
            db_config.version,
            db_config.taxon_group,
            db_config.cluster_id,
            db_config.include_singletons ?: false
        ).reference_data
    }

    emit:
    reference_data = reference_data_ch  // For silva/unite: [meta, sequences.qza, taxonomy.qza]
}

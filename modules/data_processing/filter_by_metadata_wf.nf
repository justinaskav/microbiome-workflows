
// Process to filter feature table by samples present in metadata
process FILTER_SAMPLES_BY_METADATA {
    tag "filter_samples_by_metadata"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza  // Feature table from DADA2
    path metadata_file       // Metadata file (batch1.tsv)

    output:
    path("filtered_feature_table.qza"), emit: filtered_feature_table

    script:
    """
    qiime feature-table filter-samples \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${metadata_file} \\
        --o-filtered-table filtered_feature_table.qza \\
        --verbose
    """
}

// Process to filter representative sequences based on filtered feature table
process FILTER_SEQS_BY_FILTERED_TABLE {
    tag "filter_seqs_by_filtered_table"
    conda "${params.qiime2_conda_env}"

    input:
    path rep_seqs_qza           // Representative sequences from DADA2
    path filtered_feature_table // Filtered feature table

    output:
    path("filtered_rep_seqs.qza"), emit: filtered_rep_seqs

    script:
    """
    qiime feature-table filter-seqs \\
        --i-data ${rep_seqs_qza} \\
        --i-table ${filtered_feature_table} \\
        --o-filtered-data filtered_rep_seqs.qza \\
        --verbose
    """
}

// Workflow to filter both feature table and representative sequences by metadata
workflow WF_FILTER_BY_METADATA {
    take:
    feature_table_qza    // Feature table from DADA2
    rep_seqs_qza        // Representative sequences from DADA2
    metadata_file        // Metadata file

    main:

    // Filter feature table to only include samples present in metadata
    filtered_table_ch = FILTER_SAMPLES_BY_METADATA(
        feature_table_qza,
        metadata_file
    )

    // Filter representative sequences based on filtered feature table
    filtered_seqs_ch = FILTER_SEQS_BY_FILTERED_TABLE(
        rep_seqs_qza,
        filtered_table_ch.filtered_feature_table
    )

    emit:
    filtered_feature_table = filtered_table_ch.filtered_feature_table
    filtered_rep_seqs = filtered_seqs_ch.filtered_rep_seqs
}

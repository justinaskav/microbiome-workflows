// Reusable processes for creating QIIME2 summaries and visualizations

// Process to summarize feature table
process SUMMARIZE_FEATURE_TABLE {
    tag "summarize_feature_table"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza
    path metadata_file, stageAs: "metadata.tsv"

    output:
    path("feature_table_summary.qzv"), emit: feature_table_summary_qzv

    script:
    def metadata_arg = metadata_file.name != 'NO_FILE' ? "--m-metadata-file ${metadata_file}" : ""
    """
    qiime feature-table summarize \\
        --i-table "${feature_table_qza}" \\
        --o-feature-frequencies feature_frequencies.qza \\
        --o-sample-frequencies sample_frequencies.qza \\
        --o-summary feature_table_summary.qzv \\
        ${metadata_arg} \\
        --verbose
    """
}

// Process to summarize feature table without metadata
process SUMMARIZE_FEATURE_TABLE_NO_METADATA {
    tag "summarize_feature_table_no_metadata"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza

    output:
    path("feature_table_summary.qzv"), emit: feature_table_summary_qzv

    script:
    """
    qiime feature-table summarize \\
        --i-table "${feature_table_qza}" \\
        --o-feature-frequencies feature_frequencies.qza \\
        --o-sample-frequencies sample_frequencies.qza \\
        --o-summary feature_table_summary.qzv \\
        --verbose
    """
}

// Process to tabulate representative sequences
process TABULATE_SEQS {
    tag "tabulate_rep_seqs"
    conda "${params.qiime2_conda_env}"

    input:
    path rep_seqs_qza

    output:
    path("rep_seqs_summary.qzv"), emit: rep_seqs_summary_qzv

    script:
    """
    qiime feature-table tabulate-seqs \\
        --i-data "${rep_seqs_qza}" \\
        --o-visualization rep_seqs_summary.qzv \\
        --verbose
    """
}

// Workflow for creating summaries and visualizations
workflow WF_CREATE_SUMMARIES {
    take:
        feature_table_qza
        rep_seqs_qza
        metadata_file  // Optional, can be null

    main:
        // Create feature table summary
        if (metadata_file) {
            SUMMARIZE_FEATURE_TABLE(feature_table_qza, metadata_file)
            feature_table_summary_ch = SUMMARIZE_FEATURE_TABLE.out.feature_table_summary_qzv
        } else {
            SUMMARIZE_FEATURE_TABLE_NO_METADATA(feature_table_qza)
            feature_table_summary_ch = SUMMARIZE_FEATURE_TABLE_NO_METADATA.out.feature_table_summary_qzv
        }
        
        // Create sequences summary
        TABULATE_SEQS(rep_seqs_qza)

    emit:
        feature_table_summary = feature_table_summary_ch
        rep_seqs_summary = TABULATE_SEQS.out.rep_seqs_summary_qzv
}

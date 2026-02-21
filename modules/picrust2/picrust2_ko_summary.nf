// PICRUSt2 KO Metagenome Summary Visualization
// Creates visualization of KO metagenome table

nextflow.enable.dsl = 2

process PICRUST2_KO_SUMMARY {
    tag "picrust2_ko_summary"
    conda "${params.qiime2_conda_env}"

    input:
    path ko_metagenome_qza

    output:
    path "ko_metagenome_summary.qzv", emit: ko_summary_qzv

    script:
    """
    qiime feature-table summarize \\
        --i-table "${ko_metagenome_qza}" \\
        --o-feature-frequencies feature_frequencies.qza \\
        --o-sample-frequencies sample_frequencies.qza \\
        --o-summary ko_metagenome_summary.qzv \\
        --verbose
    """
}

// Workflow wrapper for KO metagenome summary
workflow WF_PICRUST2_KO_SUMMARY {
    take:
    ko_metagenome_ch    // path to KO metagenome QZA

    main:
    summary_results = PICRUST2_KO_SUMMARY(ko_metagenome_ch)

    emit:
    ko_summary = summary_results.ko_summary_qzv
}

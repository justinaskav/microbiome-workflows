// PICRUSt2 Pathway Summary Visualization
// Creates visualization of pathway abundance table

nextflow.enable.dsl = 2

process PICRUST2_PATHWAY_SUMMARY {
    tag "picrust2_pathway_summary"
    conda "${params.qiime2_conda_env}"

    input:
    path pathway_abundance_qza

    output:
    path "pathway_abundance_summary.qzv", emit: pathway_summary_qzv

    script:
    """
    qiime feature-table summarize \\
        --i-table "${pathway_abundance_qza}" \\
        --o-feature-frequencies feature_frequencies.qza \\
        --o-sample-frequencies sample_frequencies.qza \\
        --o-summary pathway_abundance_summary.qzv \\
        --verbose
    """
}

// Workflow wrapper for pathway summary
workflow WF_PICRUST2_PATHWAY_SUMMARY {
    take:
    pathway_abundance_ch    // path to pathway abundance QZA

    main:
    summary_results = PICRUST2_PATHWAY_SUMMARY(pathway_abundance_ch)

    emit:
    pathway_summary = summary_results.pathway_summary_qzv
}

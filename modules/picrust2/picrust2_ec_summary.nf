// PICRUSt2 EC Metagenome Summary Visualization
// Creates visualization of EC metagenome table

nextflow.enable.dsl = 2

process PICRUST2_EC_SUMMARY {
    tag "picrust2_ec_summary"
    conda "${params.qiime2_conda_env}"

    input:
    path ec_metagenome_qza

    output:
    path "ec_metagenome_summary.qzv", emit: ec_summary_qzv

    script:
    """
    qiime feature-table summarize \\
        --i-table "${ec_metagenome_qza}" \\
        --o-feature-frequencies feature_frequencies.qza \\
        --o-sample-frequencies sample_frequencies.qza \\
        --o-summary ec_metagenome_summary.qzv \\
        --verbose
    """
}

// Workflow wrapper for EC metagenome summary
workflow WF_PICRUST2_EC_SUMMARY {
    take:
    ec_metagenome_ch    // path to EC metagenome QZA

    main:
    summary_results = PICRUST2_EC_SUMMARY(ec_metagenome_ch)

    emit:
    ec_summary = summary_results.ec_summary_qzv
}

// Custom taxonomy comparison workflow 
// This module provides functions to perform detailed pairwise comparison of taxonomy artifacts

// Process to compare taxonomies using a custom Python script
process COMPARE_TAXONOMIES {
    tag "compare_taxonomies"
    label 'process_low'
    
    conda "${params.qiime2_conda_env}" // Uses the QIIME2 environment

    input:
    path taxonomy1, stageAs: 'taxonomy1.qza'  // First taxonomy artifact (.qza)
    path taxonomy2, stageAs: 'taxonomy2.qza'  // Second taxonomy artifact (.qza)
    val taxonomy1_name   // Name for first taxonomy
    val taxonomy2_name   // Name for second taxonomy

    output:
    path "comparison_results.qzv", emit: comparison_results

    script:
    """
    # Run comparison script with python to ensure it works
    qiime jk visualize-taxonomy-comparison \\
        --i-tax1 taxonomy1.qza \\
        --i-tax2 taxonomy2.qza \\
        --o-visualization comparison_results.qzv \\
        --p-tax1-name "${taxonomy1_name}" \\
        --p-tax2-name "${taxonomy2_name}" \\
        --verbose
    """
}

// Workflow for comparing taxonomy artifacts
workflow COMPARE_TAXONOMY_ARTIFACTS {
    take:
    taxonomy1_ch     // Channel with first taxonomy artifact
    taxonomy2_ch     // Channel with second taxonomy artifact
    name1_ch         // Channel with name for first taxonomy
    name2_ch         // Channel with name for second taxonomy
    
    main:
    // Execute comparison
    comparison_ch = COMPARE_TAXONOMIES(
        taxonomy1_ch,
        taxonomy2_ch,
        name1_ch,
        name2_ch,
    )
    
    emit:
    comparison_results = comparison_ch.comparison_results
}

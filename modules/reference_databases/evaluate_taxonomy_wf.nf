// Unified taxonomy evaluation workflow
// Provides taxonomy evaluation functionality for both SILVA and UNITE databases

process EVALUATE_TAXONOMY {
    tag "Evaluate taxonomy - ${output_name}"
    label 'process_low'
    
    conda "${params.qiime2_conda_env}"

    input:
    path taxonomy_qzas       // One or more taxonomy artifacts (.qza)
    val output_name          // Name for the output file (without extension)
    val rank_handle_regex    // Optional regex for rank handling

    output:
    path("${output_name}.qzv"), emit: taxonomy_stats

    script:
    def taxonomies_cmd = taxonomy_qzas.collect { "--i-taxonomies ${it}" }.join(" ")
    def rank_regex_param = rank_handle_regex ? "--p-rank-handle-regex '${rank_handle_regex}'" : ""
    
    """
    qiime rescript evaluate-taxonomy \\
        ${taxonomies_cmd} \\
        ${rank_regex_param} \\
        --o-taxonomy-stats ${output_name}.qzv \\
        --verbose
    """
}

process COMPARE_TAXONOMIES {
    tag "Compare ${taxonomy1_name} vs ${taxonomy2_name}"
    label 'process_low'
    
    conda "${params.qiime2_conda_env}"

    input:
    tuple path(taxonomy1_qza), val(taxonomy1_name)
    tuple path(taxonomy2_qza), val(taxonomy2_name)
    val output_name

    output:
    path("${output_name}-comparison.qzv"), emit: comparison_stats

    script:
    """
    qiime rescript evaluate-taxonomy \\
        --i-taxonomies "${taxonomy1_qza}" "${taxonomy2_qza}" \\
        --o-taxonomy-stats "${output_name}-comparison.qzv" \\
        --verbose
    """
}

workflow EVALUATE_REFERENCE_TAXONOMY {
    take:
    taxonomy_files_ch      // Channel of taxonomy files to evaluate
    output_name           // Base name for output files
    rank_handle_regex     // Optional regex for rank handling

    main:
    // Evaluate taxonomy statistics
    evaluation_ch = EVALUATE_TAXONOMY(
        taxonomy_files_ch.collect(),
        output_name,
        rank_handle_regex ?: ''
    ).taxonomy_stats

    emit:
    evaluation = evaluation_ch
}

workflow COMPARE_REFERENCE_TAXONOMIES {
    take:
    taxonomy1_ch          // Channel: [taxonomy.qza, name]
    taxonomy2_ch          // Channel: [taxonomy.qza, name]  
    comparison_name       // Name for comparison output

    main:
    // Compare two taxonomy files
    comparison_ch = COMPARE_TAXONOMIES(
        taxonomy1_ch,
        taxonomy2_ch,
        comparison_name
    ).comparison_stats

    emit:
    comparison = comparison_ch
}

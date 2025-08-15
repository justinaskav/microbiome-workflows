// RESCRIPt evaluate-taxonomy workflow for QIIME2
// This module provides functions to evaluate taxonomic annotations using RESCRIPt

// Default parameters specific to this workflow module
// These can be overridden in nextflow.config or main script
params.rank_handle_regex = ''  // Default to empty (no regex filtering) 
params.taxonomy_evaluation_cpus = 1  // Default CPUs for the process

// Process to evaluate taxonomy artifacts using RESCRIPt
process RESCRIPT_EVALUATE_TAXONOMY {
    tag "evaluate_taxonomy"
    label 'process_low'
    
    conda "${params.qiime2_conda_env}" // Uses the QIIME2 environment with RESCRIPt installed
    cpus params.taxonomy_evaluation_cpus ?: 1

    input:
    path taxonomy_qzas       // One or more taxonomy artifacts (.qza)
    val output_name          // Name for the output file (without extension)

    output:
    path("${output_name}.qzv"), emit: taxonomy_stats_qzv

    script:
    def taxonomies_cmd = taxonomy_qzas.collect { "--i-taxonomies ${it}" }.join(" ")
    def rank_handle_regex = params.rank_handle_regex ? "--p-rank-handle-regex '${params.rank_handle_regex}'" : ""
    
    """
    qiime rescript evaluate-taxonomy \\
        ${taxonomies_cmd} \\
        ${rank_handle_regex} \\
        --o-taxonomy-stats ${output_name}.qzv \\
        --verbose
    """
}

// Workflow for evaluating taxonomies
workflow EVALUATE_TAXONOMY {
    take:
    taxonomy_qzas_ch   // Channel of one or more taxonomy artifacts (.qza)
    output_name_ch     // Channel with output name
    
    main:
    // Execute evaluation
    stats_ch = RESCRIPT_EVALUATE_TAXONOMY(
        taxonomy_qzas_ch.collect(),
        output_name_ch
    )
    
    emit:
    taxonomy_stats = stats_ch.taxonomy_stats_qzv
}

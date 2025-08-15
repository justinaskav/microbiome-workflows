// Beta diversity analysis process for R microbiome workflows

process CALC_BETA_DIVERSITY {
    tag "${tse_type}_beta_diversity"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    publishDir "${params.outputDir}/beta_diversity", mode: 'copy'
    
    input:
    tuple val(tse_type), path(tse_file), path(beta_script)
    val formula

    output:
    tuple val(tse_type), path("*_beta_results.rds"), emit: beta_results
    tuple val(tse_type), path("*_beta_plots.pdf"), emit: beta_plots, optional: true
    path "*.log", emit: log

    script:
    """
    Rscript ${beta_script} \\
        --tse ${tse_file} \\
        --tse_type ${tse_type} \\
        --formula "${formula}" \\
        --output_dir . \\
        --output_prefix ${tse_type} \\
        --methods "bray,euclidean,jaccard" \\
        --create_plots true \\
        2>&1 | tee ${tse_type}_beta_diversity.log
    """
}
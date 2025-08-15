// R Alpha Diversity Module
nextflow.enable.dsl = 2

process CALC_ALPHA_DIVERSITY {
    tag "calc_alpha_diversity_${tse_type}"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    label 'process_medium'

    input:
    tuple val(tse_type), path(tse_file), path(calc_alpha_script)
    val alpha_indices
    val alpha_names  
    val niter

    output:
    path "tse_${tse_type}_with_alpha.rds", emit: tse_with_alpha
    path "alpha_diversity_${tse_type}.rds", emit: alpha_results

    script:
    """
    Rscript ${calc_alpha_script} \
        --tse ${tse_file} \
        --tse_type ${tse_type} \
        --output_dir . \
        --indices "${alpha_indices}" \
        --names "${alpha_names}" \
        --niter ${niter}
    """
}

workflow ALPHA_DIVERSITY_WF {
    take:
    tse_data_with_script   // Channel with [type, tse_file, script]
    alpha_indices         // Diversity indices to calculate
    alpha_names          // Names for indices
    niter               // Number of iterations

    main:
    alpha_results = CALC_ALPHA_DIVERSITY(
        tse_data_with_script,
        alpha_indices,
        alpha_names,
        niter
    )

    emit:
    tse_with_alpha = alpha_results.tse_with_alpha
    alpha_results = alpha_results.alpha_results
}
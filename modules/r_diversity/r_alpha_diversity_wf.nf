// R-based alpha diversity calculation for TSE objects

process CALC_R_ALPHA_DIVERSITY {
    tag "calc_alpha_diversity_${tse_type}"
    conda params.r_conda_env

    input:
    tuple val(tse_type), path(tse_file), val(indices), val(names), val(niter), val(seed), path(calc_alpha_script), path(utils_script)

    output:
    path "tse_${tse_type}_with_alpha.rds", emit: tse_with_alpha
    path "alpha_diversity_${tse_type}.rds", emit: alpha_results

    script:
    """
    Rscript ${calc_alpha_script} \
        --tse ${tse_file} \
        --tse_type ${tse_type} \
        --output_dir . \
        --indices ${indices} \
        --names ${names} \
        --niter ${niter} \
        --seed ${seed}
    """
}

workflow WF_R_ALPHA_DIVERSITY {
    take:
    tse_data   // Channel of [tse_type, tse_file]
    indices    // String: comma-separated alpha diversity indices
    names      // String: comma-separated display names
    niter      // Integer: number of rarefaction iterations
    seed       // Integer: random seed

    main:
    def calc_alpha_script = file("${moduleDir}/../../scripts/calc_alpha_diversity.R")
    def utils_script = file("${moduleDir}/../../scripts/r_utils.R")

    process_input = tse_data.map { tse_type, tse_file ->
        [tse_type, tse_file, indices, names, niter, seed, calc_alpha_script, utils_script]
    }

    CALC_R_ALPHA_DIVERSITY(process_input)

    emit:
    tse_with_alpha = CALC_R_ALPHA_DIVERSITY.out.tse_with_alpha
    alpha_results = CALC_R_ALPHA_DIVERSITY.out.alpha_results
}

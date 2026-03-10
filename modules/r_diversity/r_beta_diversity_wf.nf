// R-based beta diversity calculation for TSE objects

process CALC_R_BETA_DIVERSITY {
    tag "calc_beta_diversity_${tse_type}"
    conda params.r_conda_env

    input:
    tuple val(tse_type), path(tse_file), val(distance_methods), val(ncomponents), val(niter), val(seed), path(calc_beta_script), path(utils_script)

    output:
    path "tse_${tse_type}_with_beta.rds", emit: tse_with_beta

    script:
    """
    Rscript ${calc_beta_script} \
        --tse ${tse_file} \
        --tse_type ${tse_type} \
        --output_dir . \
        --distance_methods ${distance_methods} \
        --ncomponents ${ncomponents} \
        --niter ${niter} \
        --seed ${seed}
    """
}

workflow WF_R_BETA_DIVERSITY {
    take:
    tse_data          // Channel of [tse_type, tse_file]
    distance_methods  // String: comma-separated distance methods
    ncomponents       // Integer: number of MDS components
    niter             // Integer: number of rarefaction iterations
    seed              // Integer: random seed

    main:
    def calc_beta_script = file("${moduleDir}/../../scripts/calc_beta_diversity.R")
    def utils_script = file("${moduleDir}/../../scripts/r_utils.R")

    process_input = tse_data.map { tse_type, tse_file ->
        [tse_type, tse_file, distance_methods, ncomponents, niter, seed, calc_beta_script, utils_script]
    }

    CALC_R_BETA_DIVERSITY(process_input)

    emit:
    tse_with_beta = CALC_R_BETA_DIVERSITY.out.tse_with_beta
}

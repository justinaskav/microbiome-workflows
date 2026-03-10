// R Beta Diversity Subworkflow
// Calculates beta diversity distance matrices and MDS ordination on TSE objects

nextflow.enable.dsl = 2

include { WF_R_BETA_DIVERSITY } from '../modules/r_diversity/r_beta_diversity_wf.nf'

workflow R_BETA_DIVERSITY {
    take:
    tse_data      // Channel of [tse_type, tse_file]
    beta_options  // Map with distance_methods, ncomponents (can be empty/null)

    main:
    def config = [
        distance_methods: beta_options?.distance_methods ?: params.r_diversity?.beta?.distance_methods ?: 'jaccard,bray,unifrac_unweighted,unifrac_weighted',
        ncomponents:      beta_options?.ncomponents      ?: params.r_diversity?.beta?.ncomponents      ?: 2,
        niter:            beta_options?.niter             ?: params.r_diversity?.niter                  ?: 500,
        seed:             beta_options?.seed              ?: params.r_diversity?.seed                   ?: 14100715,
    ]

    WF_R_BETA_DIVERSITY(
        tse_data,
        config.distance_methods,
        config.ncomponents,
        config.niter,
        config.seed,
    )

    emit:
    tse_with_beta = WF_R_BETA_DIVERSITY.out.tse_with_beta
}

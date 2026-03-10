// R Alpha Diversity Subworkflow
// Calculates alpha diversity indices on TSE objects using mia::addAlpha

nextflow.enable.dsl = 2

include { WF_R_ALPHA_DIVERSITY } from '../modules/r_diversity/r_alpha_diversity_wf.nf'

workflow R_ALPHA_DIVERSITY {
    take:
    tse_data       // Channel of [tse_type, tse_file]
    alpha_options  // Map with indices, names (can be empty/null)

    main:
    def config = [
        indices: alpha_options?.indices ?: params.r_diversity?.alpha?.indices ?: 'shannon,faith_diversity',
        names:   alpha_options?.names   ?: params.r_diversity?.alpha?.names   ?: 'Shannon,Faith_PD',
        niter:   alpha_options?.niter   ?: params.r_diversity?.niter          ?: 500,
        seed:    alpha_options?.seed    ?: params.r_diversity?.seed           ?: 14100715,
    ]

    WF_R_ALPHA_DIVERSITY(
        tse_data,
        config.indices,
        config.names,
        config.niter,
        config.seed,
    )

    emit:
    tse_with_alpha = WF_R_ALPHA_DIVERSITY.out.tse_with_alpha
    alpha_results = WF_R_ALPHA_DIVERSITY.out.alpha_results
}

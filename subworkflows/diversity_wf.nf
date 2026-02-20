// Diversity Subworkflow
// Runs diversity core metrics analysis (phylogenetic or non-phylogenetic based on amplicon type)
// Supports explicit sampling depth or auto-inference from the feature table.

nextflow.enable.dsl = 2

include { WF_DIVERSITY_CORE_METRICS ; WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC } from '../modules/diversity/diversity_core_metrics_wf.nf'
include { WF_INFER_SAMPLING_DEPTH } from '../modules/diversity/infer_sampling_depth.nf'

workflow DIVERSITY {
    take:
    feature_table      // Filtered feature table
    rooted_tree        // Rooted phylogenetic tree (can be empty for ITS)
    metadata_file      // Metadata file
    diversity_options  // Map with sampling_depth, amplicon_type, database, skip, auto_depth_percentile

    main:

    // Initialize empty channels for all outputs
    rarefied_table_ch = Channel.empty()
    jaccard_emperor_ch = Channel.empty()
    bray_curtis_emperor_ch = Channel.empty()

    if (!diversity_options?.skip) {
        // Step 1: Resolve sampling depth
        if (diversity_options?.sampling_depth) {
            sampling_depth_ch = Channel.value(diversity_options.sampling_depth)
        } else {
            def percentile = diversity_options?.auto_depth_percentile ?: 0.0
            inferred = WF_INFER_SAMPLING_DEPTH(feature_table, percentile, metadata_file)
            sampling_depth_ch = inferred.sampling_depth_file
                .map { it.text.trim().toInteger() }
                .first()
        }

        // Step 2: Run appropriate diversity method
        if (diversity_options?.amplicon_type == "ITS" || diversity_options?.database == "unite") {
            diversity_results = WF_DIVERSITY_CORE_METRICS(
                feature_table,
                metadata_file,
                sampling_depth_ch,
                [id: "diversity_analysis"],
            )
        } else {
            diversity_results = WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC(
                feature_table,
                rooted_tree,
                metadata_file,
                sampling_depth_ch,
                [id: "diversity_analysis"],
            )
        }

        rarefied_table_ch = diversity_results.rarefied_table
        jaccard_emperor_ch = diversity_results.jaccard_emperor
        bray_curtis_emperor_ch = diversity_results.bray_curtis_emperor
    }

    emit:
    rarefied_table = rarefied_table_ch
    jaccard_emperor = jaccard_emperor_ch
    bray_curtis_emperor = bray_curtis_emperor_ch
}

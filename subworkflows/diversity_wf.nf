// Diversity Subworkflow
// Runs diversity core metrics analysis (phylogenetic or non-phylogenetic based on amplicon type)

nextflow.enable.dsl = 2

include { WF_DIVERSITY_CORE_METRICS ; WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC } from '../modules/diversity/diversity_core_metrics_wf.nf'

workflow DIVERSITY {
    take:
    feature_table      // Filtered feature table
    rooted_tree        // Rooted phylogenetic tree (can be empty for ITS)
    metadata_file      // Metadata file
    diversity_options  // Map with sampling_depth, amplicon_type, database

    main:

    // Initialize empty channels for all outputs
    rarefied_table_ch = Channel.empty()
    jaccard_emperor_ch = Channel.empty()
    bray_curtis_emperor_ch = Channel.empty()

    if (diversity_options?.sampling_depth) {
        // Use non-phylogenetic diversity for ITS (fungal) data, phylogenetic for 16S (bacterial) data
        if (diversity_options?.amplicon_type == "ITS" || diversity_options?.database == "unite") {
            diversity_results = WF_DIVERSITY_CORE_METRICS(
                feature_table,
                metadata_file,
                diversity_options.sampling_depth,
                [id: "diversity_analysis"],
            )
            rarefied_table_ch = diversity_results.rarefied_table
            jaccard_emperor_ch = diversity_results.jaccard_emperor
            bray_curtis_emperor_ch = diversity_results.bray_curtis_emperor
        }
        else {
            // Default phylogenetic diversity for 16S/bacterial data
            diversity_results = WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC(
                feature_table,
                rooted_tree,
                metadata_file,
                diversity_options.sampling_depth,
                [id: "diversity_analysis"],
            )
            rarefied_table_ch = diversity_results.rarefied_table
            jaccard_emperor_ch = diversity_results.jaccard_emperor
            bray_curtis_emperor_ch = diversity_results.bray_curtis_emperor
        }
    }

    emit:
    rarefied_table = rarefied_table_ch
    jaccard_emperor = jaccard_emperor_ch
    bray_curtis_emperor = bray_curtis_emperor_ch
}

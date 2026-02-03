// PICRUSt2 Subworkflow
// Orchestrates PICRUSt2 functional prediction, visualization, and diversity analysis
// Supports explicit pathway diversity sampling depth or auto-inference from pathway abundance.

nextflow.enable.dsl = 2

include { WF_PICRUST2_FULL_PIPELINE } from '../modules/picrust2/picrust2_full_pipeline.nf'
include { WF_PICRUST2_KO_SUMMARY } from '../modules/picrust2/picrust2_ko_summary.nf'
include { WF_PICRUST2_EC_SUMMARY } from '../modules/picrust2/picrust2_ec_summary.nf'
include { WF_PICRUST2_PATHWAY_SUMMARY } from '../modules/picrust2/picrust2_pathway_summary.nf'
include { WF_PICRUST2_PATHWAY_DIVERSITY } from '../modules/picrust2/picrust2_pathway_diversity.nf'
include { WF_INFER_SAMPLING_DEPTH as WF_INFER_PATHWAY_DEPTH } from '../modules/diversity/infer_sampling_depth.nf'

workflow PICRUST2 {
    take:
    feature_table      // Filtered feature table from FILTERING step
    rep_seqs           // Filtered representative sequences from FILTERING step
    metadata_file      // Metadata file
    picrust2_options   // Map with threads, diversity_sampling_depth, skip_pathway_diversity, auto_depth_percentile

    main:

    // Configure PICRUSt2 with fallback to DADA2 threads and global defaults
    picrust2_config = [
        threads: picrust2_options?.threads ?: params.dada2?.n_threads ?: 8,
        diversity_sampling_depth: picrust2_options?.diversity_sampling_depth ?: params.picrust2?.diversity_sampling_depth,
        skip_pathway_diversity: picrust2_options?.skip_pathway_diversity ?: false,
        auto_depth_percentile: picrust2_options?.auto_depth_percentile ?: 0.0,
    ]

    // Step 1: Run PICRUSt2 full pipeline
    picrust2_results = WF_PICRUST2_FULL_PIPELINE(
        feature_table,
        rep_seqs,
        picrust2_config.threads
    )

    // Step 2: Create summary visualizations
    ko_summary_results = WF_PICRUST2_KO_SUMMARY(
        picrust2_results.ko_metagenome
    )

    ec_summary_results = WF_PICRUST2_EC_SUMMARY(
        picrust2_results.ec_metagenome
    )

    pathway_summary_results = WF_PICRUST2_PATHWAY_SUMMARY(
        picrust2_results.pathway_abundance
    )

    // Step 3: Conditional diversity analysis on pathway abundances
    // Initialize empty channels
    pathway_rarefied_table_ch = Channel.empty()
    pathway_jaccard_emperor_ch = Channel.empty()
    pathway_bray_curtis_emperor_ch = Channel.empty()
    pathway_observed_features_ch = Channel.empty()
    pathway_shannon_ch = Channel.empty()
    pathway_evenness_ch = Channel.empty()
    pathway_jaccard_dm_ch = Channel.empty()
    pathway_bray_curtis_dm_ch = Channel.empty()
    pathway_jaccard_pcoa_ch = Channel.empty()
    pathway_bray_curtis_pcoa_ch = Channel.empty()

    if (!picrust2_config.skip_pathway_diversity) {
        // Resolve sampling depth
        if (picrust2_config.diversity_sampling_depth) {
            sampling_depth_ch = Channel.value(picrust2_config.diversity_sampling_depth)
        } else {
            inferred = WF_INFER_PATHWAY_DEPTH(
                picrust2_results.pathway_abundance,
                picrust2_config.auto_depth_percentile
            )
            sampling_depth_ch = inferred.sampling_depth_file
                .map { it.text.trim().toInteger() }
                .first()
        }

        diversity_results = WF_PICRUST2_PATHWAY_DIVERSITY(
            picrust2_results.pathway_abundance,
            metadata_file,
            sampling_depth_ch
        )

        pathway_rarefied_table_ch = diversity_results.rarefied_table
        pathway_jaccard_emperor_ch = diversity_results.jaccard_emperor
        pathway_bray_curtis_emperor_ch = diversity_results.bray_curtis_emperor
        pathway_observed_features_ch = diversity_results.observed_features_vector
        pathway_shannon_ch = diversity_results.shannon_vector
        pathway_evenness_ch = diversity_results.evenness_vector
        pathway_jaccard_dm_ch = diversity_results.jaccard_distance_matrix
        pathway_bray_curtis_dm_ch = diversity_results.bray_curtis_distance_matrix
        pathway_jaccard_pcoa_ch = diversity_results.jaccard_pcoa_results
        pathway_bray_curtis_pcoa_ch = diversity_results.bray_curtis_pcoa_results
    }

    emit:
    // PICRUSt2 functional predictions
    ko_metagenome = picrust2_results.ko_metagenome
    ec_metagenome = picrust2_results.ec_metagenome
    pathway_abundance = picrust2_results.pathway_abundance

    // Summary visualizations
    ko_summary = ko_summary_results.ko_summary
    ec_summary = ec_summary_results.ec_summary
    pathway_summary = pathway_summary_results.pathway_summary

    // Pathway diversity results (may be empty if explicitly skipped)
    pathway_rarefied_table = pathway_rarefied_table_ch
    pathway_jaccard_emperor = pathway_jaccard_emperor_ch
    pathway_bray_curtis_emperor = pathway_bray_curtis_emperor_ch
    pathway_observed_features = pathway_observed_features_ch
    pathway_shannon = pathway_shannon_ch
    pathway_evenness = pathway_evenness_ch
    pathway_jaccard_dm = pathway_jaccard_dm_ch
    pathway_bray_curtis_dm = pathway_bray_curtis_dm_ch
    pathway_jaccard_pcoa = pathway_jaccard_pcoa_ch
    pathway_bray_curtis_pcoa = pathway_bray_curtis_pcoa_ch
}

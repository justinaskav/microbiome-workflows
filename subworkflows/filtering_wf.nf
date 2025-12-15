// Filtering Subworkflow
// Performs metadata filtering and optional decontamination (frequency or prevalence)

nextflow.enable.dsl = 2

include { WF_FILTER_BY_METADATA } from '../modules/data_processing/filter_by_metadata_wf.nf'
include { WF_DECONTAM_FREQUENCY ; WF_DECONTAM_PREVALENCE } from '../modules/decontam/decontam_wf.nf'

workflow FILTERING {
    take:
    feature_table      // Feature table from DADA2
    rep_seqs           // Representative sequences from DADA2
    metadata_file      // Metadata file
    decontam_options   // Map with method, freq_concentration_column, prev_control_column, prev_control_indicator, threshold

    main:

    decontam_config = [
        method: decontam_options?.method ?: params.decontam?.method ?: "none",
        freq_concentration_column: decontam_options?.freq_concentration_column ?: params.decontam?.freq_concentration_column,
        prev_control_column: decontam_options?.prev_control_column ?: params.decontam?.prev_control_column,
        prev_control_indicator: decontam_options?.prev_control_indicator ?: params.decontam?.prev_control_indicator,
        threshold: decontam_options?.threshold ?: params.decontam?.threshold ?: 0.1,
    ]

    // Filter DADA2 outputs to only include samples present in metadata
    filtering_results_ch = WF_FILTER_BY_METADATA(
        feature_table,
        rep_seqs,
        metadata_file,
    )

    // Conditional decontamination based on method specified
    if (decontam_config.method == "frequency") {
        decontam_results_ch = WF_DECONTAM_FREQUENCY(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_file,
            decontam_config.freq_concentration_column,
            decontam_config.threshold,
        )

        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores = decontam_results_ch.decontam_scores
        decontam_scores_viz = decontam_results_ch.decontam_scores_viz
    }
    else if (decontam_config.method == "prevalence") {
        decontam_results_ch = WF_DECONTAM_PREVALENCE(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_file,
            decontam_config.prev_control_column,
            decontam_config.prev_control_indicator,
            decontam_config.threshold,
        )

        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores = decontam_results_ch.decontam_scores
        decontam_scores_viz = decontam_results_ch.decontam_scores_viz
    }
    else {
        // No decontamination - use metadata-filtered results directly
        filtered_feature_table = filtering_results_ch.filtered_feature_table
        filtered_rep_seqs = filtering_results_ch.filtered_rep_seqs
        decontam_scores = Channel.empty()
        decontam_scores_viz = Channel.empty()
    }

    emit:
    filtered_feature_table = filtered_feature_table
    filtered_rep_seqs = filtered_rep_seqs
    decontam_scores = decontam_scores
    decontam_scores_viz = decontam_scores_viz
}

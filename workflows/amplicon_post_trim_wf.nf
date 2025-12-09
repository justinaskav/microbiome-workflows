// Amplicon Post-Trim Workflow
// Takes trimmed reads (potentially merged from multiple batches) and runs DADA2 → downstream

nextflow.enable.dsl = 2

include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'
include { WF_FILTER_BY_METADATA } from '../modules/data_processing/filter_by_metadata_wf.nf'
include { WF_DECONTAM_FREQUENCY ; WF_DECONTAM_PREVALENCE } from '../modules/decontam/decontam_wf.nf'
include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'
include { EVALUATE_REFERENCE_TAXONOMY } from '../modules/reference_databases/evaluate_taxonomy_wf.nf'
include { COMPARE_TAXONOMY_ARTIFACTS } from '../modules/taxonomy/taxonomy_comparison_wf.nf'
include { WF_EVALUATE_GRAM_STAINING } from '../modules/taxonomy/evaluate_gram_staining_wf.nf'
include { WF_DIVERSITY_CORE_METRICS ; WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC } from '../modules/diversity/diversity_core_metrics_wf.nf'

workflow AMPLICON_POST_TRIM {
    take:
    trimmed_reads_ch      // Channel of [sample_id, [r1, r2]] - already trimmed
    metadata_file         // Merged metadata file
    dada2_options
    decontam_options
    taxonomy_options
    database_ch_list
    other_options

    main:

    dada2_config = [
        trunc_len_f: dada2_options.trunc_len_f,
        trunc_len_r: dada2_options.trunc_len_r,
        trim_left_f: dada2_options.trim_left_f ?: 0,
        trim_left_r: dada2_options.trim_left_r ?: 0,
        n_reads_learn: dada2_options.n_reads_learn ?: params.dada2?.n_reads_learn ?: 1000000,
        n_threads: dada2_options.n_threads ?: params.dada2?.n_threads ?: 8,
        sample_prefix_to_remove: dada2_options.sample_prefix_to_remove ?: params.dada2?.sample_prefix_to_remove,
    ]

    decontam_config = [
        method: decontam_options?.method ?: params.decontam?.method ?: "none",
        freq_concentration_column: decontam_options?.freq_concentration_column ?: params.decontam?.freq_concentration_column,
        prev_control_column: decontam_options?.prev_control_column ?: params.decontam?.prev_control_column,
        prev_control_indicator: decontam_options?.prev_control_indicator ?: params.decontam?.prev_control_indicator,
        threshold: decontam_options?.threshold ?: params.decontam?.threshold ?: 0.1,
    ]

    taxonomy_config = [
        enable_dada2_taxonomy: taxonomy_options?.enable_dada2_taxonomy ?: params.taxonomy?.enable_dada2_taxonomy ?: false,
        enable_comparison: taxonomy_options?.enable_comparison ?: params.taxonomy?.enable_comparison ?: false,
    ]

    // DADA2 denoising on merged/trimmed reads
    dada2_output_ch = WF_DENOISE_DADA2(
        trimmed_reads_ch,
        dada2_config.trunc_len_f,
        dada2_config.trunc_len_r,
        dada2_config.trim_left_f,
        dada2_config.trim_left_r,
        dada2_config.n_reads_learn,
        dada2_config.n_threads,
        dada2_config.sample_prefix_to_remove,
    )

    dada2_feature_table = dada2_output_ch.feature_table
    dada2_rep_seqs = dada2_output_ch.representative_sequences

    // Filter to samples in metadata
    filtering_results_ch = WF_FILTER_BY_METADATA(
        dada2_feature_table,
        dada2_rep_seqs,
        metadata_file,
    )

    // Decontamination
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
    } else if (decontam_config.method == "prevalence") {
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
    } else {
        filtered_feature_table = filtering_results_ch.filtered_feature_table
        filtered_rep_seqs = filtering_results_ch.filtered_rep_seqs
        decontam_scores = Channel.empty()
        decontam_scores_viz = Channel.empty()
    }

    // Taxonomy
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        filtered_rep_seqs,
        database_ch_list.classifier,
        filtered_feature_table,
        metadata_file,
    )

    classified_taxa_ch = taxonomy_results_ch.classification

    // Optional DADA2 taxonomy
    if (taxonomy_config.enable_dada2_taxonomy) {
        dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
            filtered_rep_seqs,
            filtered_feature_table,
            database_ch_list.refs.trainset,
            database_ch_list.refs.species,
            metadata_file
        )

        if (taxonomy_config.enable_comparison) {
            taxonomy_eval_ch = EVALUATE_REFERENCE_TAXONOMY(
                classified_taxa_ch.mix(dada2_taxonomy_results_ch.qiime2_taxonomy),
                "taxonomy_evaluation",
                "",
            ).evaluation

            taxonomy_comparison_ch = COMPARE_TAXONOMY_ARTIFACTS(
                classified_taxa_ch,
                dada2_taxonomy_results_ch.qiime2_taxonomy,
                Channel.value("Feature-Classifier"),
                Channel.value("DADA2"),
            ).comparison_results

            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
        } else {
            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
            taxonomy_eval_ch = Channel.empty()
            taxonomy_comparison_ch = Channel.empty()
        }
    } else {
        dada2_taxonomy_ch = Channel.empty()
        taxonomy_eval_ch = Channel.empty()
        taxonomy_comparison_ch = Channel.empty()
    }

    // Gram staining
    WF_EVALUATE_GRAM_STAINING(classified_taxa_ch, filtered_feature_table, metadata_file)

    // Phylogeny
    ALIGN_AND_TREE_QIIME2(filtered_rep_seqs)

    // Diversity (optional)
    if (other_options?.diversity_sampling_depth) {
        if (other_options?.amplicon_type == "ITS" || other_options?.database == "unite") {
            diversity_results_ch = WF_DIVERSITY_CORE_METRICS(
                filtered_feature_table,
                metadata_file,
                other_options.diversity_sampling_depth,
                [id: "diversity_analysis"],
            )
        } else {
            diversity_results_ch = WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC(
                filtered_feature_table,
                ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza,
                metadata_file,
                other_options.diversity_sampling_depth,
                [id: "diversity_analysis"],
            )
        }
    } else {
        diversity_results_ch = Channel.empty()
    }

    emit:
    feature_table = filtered_feature_table
    rep_seqs = filtered_rep_seqs
    taxonomy = classified_taxa_ch
    taxonomy_summary = taxonomy_results_ch.classification_summary
    dada2_taxonomy = dada2_taxonomy_ch
    taxonomy_evaluation = taxonomy_eval_ch
    taxonomy_comparison = taxonomy_comparison_ch
    decontam_scores = decontam_scores
    decontam_scores_viz = decontam_scores_viz
    phylogeny_rooted = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
    phylogeny_unrooted = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
    gram_staining_evaluation = WF_EVALUATE_GRAM_STAINING.out.gram_staining_eval_results
    diversity_results = diversity_results_ch
    dada2_ft_summary = dada2_output_ch.feature_table_summary
    dada2_rep_seqs_summary = dada2_output_ch.rep_seqs_summary
    dada2_denoising_stats_summary = dada2_output_ch.denoising_stats_summary
    dada2_base_transition_stats = dada2_output_ch.base_transition_stats
}

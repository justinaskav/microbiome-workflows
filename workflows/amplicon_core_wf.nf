// Core Amplicon Processing Workflow  
// Clean workflow with no input format branching - takes standardized inputs

nextflow.enable.dsl = 2

// Include required workflows and modules
include { WF_QUALITY_CONTROL ; WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_CUTADAPT } from '../modules/quality_control/quality_control_wf.nf'
include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'
include { WF_FILTER_BY_METADATA } from '../modules/data_processing/filter_by_metadata_wf.nf'
include { WF_DECONTAM_FREQUENCY ; WF_DECONTAM_PREVALENCE } from '../modules/decontam/decontam_wf.nf'
include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { DOWNLOAD_REFERENCE_DATABASE } from '../modules/reference_databases/download_reference_db_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'
include { EVALUATE_REFERENCE_TAXONOMY } from '../modules/reference_databases/evaluate_taxonomy_wf.nf'
include { COMPARE_TAXONOMY_ARTIFACTS } from '../modules/taxonomy/taxonomy_comparison_wf.nf'
include { WF_EVALUATE_GRAM_STAINING } from '../modules/taxonomy/evaluate_gram_staining_wf.nf'
include { WF_DIVERSITY_CORE_METRICS ; WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC } from '../modules/diversity/diversity_core_metrics_wf.nf'

// Core amplicon processing workflow with standardized inputs
workflow AMPLICON_CORE {
    take:
    read_pairs_ch // Standardized read pairs channel
    metadata_path // Metadata 
    barcode_primer_info_ch // Barcode/primer info (may be empty)
    primers
    trimming_options // Unified cutadapt configuration
    dada2_options // Unified DADA2 configuration
    decontam_options // Unified decontamination configuration
    taxonomy_options // Unified taxonomy configuration
    database_ch_list // Classifier and reference databases (SILVA or UNITE)
    other_options // Other options (e.g., sample prefix to remove, diversity settings)

    main: 

    dada2_config = [
        trunc_len_f: dada2_options.trunc_len_f,
        trunc_len_r: dada2_options.trunc_len_r,
        trim_left_f: dada2_options.trim_left_f,
        trim_left_r: dada2_options.trim_left_r,
        n_reads_learn: dada2_options.n_reads_learn ?: params.dada2?.n_reads_learn ?: 1000000,
        n_threads: dada2_options.n_threads ?: params.dada2?.n_threads ?: 8,
        sample_prefix_to_remove: dada2_options.sample_prefix_to_remove ?: params.dada2?.sample_prefix_to_remove,
    ]

    decontam_config = [
        method: decontam_options?.method ?: params.decontam?.method ?: "none",
        freq_concentration_column: decontam_options?.freq_concentration_column ?: params.decontam?.freq_concentration_column ?: null,
        prev_control_column: decontam_options?.prev_control_column ?: params.decontam?.prev_control_column ?: null,
        prev_control_indicator: decontam_options?.prev_control_indicator ?: params.decontam?.prev_control_indicator ?: null,
        threshold: decontam_options?.threshold ?: params.decontam?.threshold ?: 0.1,
    ]

    taxonomy_config = [
        enable_dada2_taxonomy: taxonomy_options?.enable_dada2_taxonomy ?: params.taxonomy?.enable_dada2_taxonomy ?: false,
        enable_comparison: taxonomy_options?.enable_comparison ?: params.taxonomy?.enable_comparison ?: false,
    ]

    trimming_config = [
        enabled: trimming_options.enabled ?: true,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
        skip_polyg_trim: trimming_options.skip_polyg_trim ?: false,
    ]

    // Run regular quality control on original reads for comparison
    WF_QUALITY_CONTROL(read_pairs_ch)

    trimming_results = WF_CUTADAPT(
        read_pairs_ch,
        barcode_primer_info_ch,
        trimming_config,
        primers
    )

    // Process exported fastqs for QC after cutadapt
    trimming_results.exported_fastqs
        .map { sample_name, fastq1, fastq2 ->
            return tuple(sample_name, [fastq1, fastq2])
        }
        .set { exported_fastqs_pairs_ch }

    // Use exported fastqs for downstream analysis (DADA2 needs raw FASTQ files)
    trimmed_sequences_for_dada2 = exported_fastqs_pairs_ch

    // Run QC on the exported fastqs after cutadapt
    WF_QUALITY_CONTROL_AFTER_CUTADAPT(exported_fastqs_pairs_ch)

    // DADA2 denoising
    dada2_output_ch = WF_DENOISE_DADA2(
        trimmed_sequences_for_dada2,
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
    dada2_denoising_stats = dada2_output_ch.denoising_stats

    // Filter DADA2 outputs to only include samples present in metadata
    filtering_results_ch = WF_FILTER_BY_METADATA(
        dada2_feature_table,
        dada2_rep_seqs,
        metadata_path,
    )

    // Conditional decontamination based on method specified
    if (decontam_config.method == "frequency") {
        decontam_results_ch = WF_DECONTAM_FREQUENCY(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_path,
            decontam_config.freq_concentration_column,
            decontam_config.threshold,
        )

        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores_freq = decontam_results_ch.decontam_scores
        decontam_scores_freq_viz = decontam_results_ch.decontam_scores_viz
        decontam_scores_prev = Channel.empty()
        decontam_scores_prev_viz = Channel.empty()
    }
    else if (decontam_config.method == "prevalence") {
        decontam_results_ch = WF_DECONTAM_PREVALENCE(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_path,
            decontam_config.prev_control_column,
            decontam_config.prev_control_indicator,
            decontam_config.threshold,
        )

        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores_freq = Channel.empty()
        decontam_scores_freq_viz = Channel.empty()
        decontam_scores_prev = decontam_results_ch.decontam_scores
        decontam_scores_prev_viz = decontam_results_ch.decontam_scores_viz
    }
    else {
        // No decontamination - use metadata-filtered results directly
        filtered_feature_table = filtering_results_ch.filtered_feature_table
        filtered_rep_seqs = filtering_results_ch.filtered_rep_seqs
        decontam_scores_freq = Channel.empty()
        decontam_scores_freq_viz = Channel.empty()
        decontam_scores_prev = Channel.empty()
        decontam_scores_prev_viz = Channel.empty()
    }

    // Taxonomic Classification with feature-classifier 
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        filtered_rep_seqs,
        database_ch_list.classifier,
        filtered_feature_table,
        metadata_path,
    )

    classified_taxa_ch = taxonomy_results_ch.classification
    classified_taxa_summary_ch = taxonomy_results_ch.classification_summary

    if (taxonomy_config.enable_dada2_taxonomy) {
        dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
            filtered_rep_seqs,
            filtered_feature_table,
            database_ch_list.refs.trainset,
            database_ch_list.refs.species,
        )

        if (taxonomy_config.enable_comparison) {

            // Evaluate taxonomic classifications using RESCRIPt
            taxonomy_eval_ch = EVALUATE_REFERENCE_TAXONOMY(
                classified_taxa_ch.mix(dada2_taxonomy_results_ch.qiime2_taxonomy),
                "taxonomy_evaluation",
                "",
            ).evaluation

            taxonomy_comparison_ch_results = COMPARE_TAXONOMY_ARTIFACTS(
                classified_taxa_ch,
                dada2_taxonomy_results_ch.qiime2_taxonomy,
                Channel.value("Feature-Classifier"),
                Channel.value("DADA2"),
            )

            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
            taxonomy_evaluation_ch = taxonomy_eval_ch
            taxonomy_comparison_ch = taxonomy_comparison_ch_results.comparison_results
        }
        else {
            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
            taxonomy_evaluation_ch = Channel.empty()
            taxonomy_comparison_ch = Channel.empty()
        }
    }
    else {
        // Skip DADA2 taxonomy - create empty channels
        dada2_taxonomy_ch = Channel.empty()
        taxonomy_evaluation_ch = Channel.empty()
        taxonomy_comparison_ch = Channel.empty()
    }

    // Gram Staining Evaluation (Optional)
    WF_EVALUATE_GRAM_STAINING(
        classified_taxa_ch,
        filtered_feature_table,
        metadata_path,
    )

    // Phylogenetic Analysis (Optional)
    ALIGN_AND_TREE_QIIME2(filtered_rep_seqs)

    // Diversity Analysis (Optional)
    diversity_results_ch = Channel.empty()
    if (other_options?.diversity_sampling_depth) {
        // Use phylogenetic diversity for 16S (bacterial) data, non-phylogenetic for ITS (fungal) data
        if (other_options?.amplicon_type == "ITS" || other_options?.database == "unite") {
            diversity_results_ch = WF_DIVERSITY_CORE_METRICS(
                filtered_feature_table,
                metadata_path,
                other_options.diversity_sampling_depth,
                [id: "diversity_analysis"],
            )
        }
        else {
            // Default phylogenetic diversity for 16S/bacterial data
            diversity_results_ch = WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC(
                filtered_feature_table,
                ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza,
                metadata_path,
                other_options.diversity_sampling_depth,
                [id: "diversity_analysis"],
            )
        }
    }

    emit:
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    fastqc_results_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.fastqc_results
    multiqc_report_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.multiqc_report
    dada2_denoising_stats = dada2_denoising_stats
    dada2_ft_summary = WF_DENOISE_DADA2.out.feature_table_summary
    dada2_rep_seqs_summary = WF_DENOISE_DADA2.out.rep_seqs_summary
    dada2_denoising_stats_summary = WF_DENOISE_DADA2.out.denoising_stats_summary
    decontam_scores_freq = decontam_scores_freq
    decontam_scores_freq_viz = decontam_scores_freq_viz
    decontam_scores_prev = decontam_scores_prev
    decontam_scores_prev_viz = decontam_scores_prev_viz
    taxonomy = classified_taxa_ch
    taxonomy_summary = classified_taxa_summary_ch
    phylogeny_unrooted = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
    phylogeny_rooted = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
    dada2_taxonomy = dada2_taxonomy_ch
    taxonomy_evaluation = taxonomy_evaluation_ch
    taxonomy_comparison = taxonomy_comparison_ch
    gram_staining_evaluation = WF_EVALUATE_GRAM_STAINING.out.gram_staining_eval_results
    diversity_results = diversity_results_ch
    feature_table = filtered_feature_table
    rep_seqs = filtered_rep_seqs
}

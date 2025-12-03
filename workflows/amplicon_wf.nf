// Include input processing workflow
include { WF_INPUT_PROCESSING } from './amplicon_input_processing_wf.nf'

// Include core processing workflows
include { AMPLICON_CORE } from './amplicon_core_wf.nf'
include { AMPLICON_QC_ONLY } from './amplicon_qc_only_wf.nf'

// Multi-modal workflow dispatcher
workflow AMPLICON_WORKFLOW {
    take:
    name
    metadata_path
    input_options
    primers
    trimming_options
    dada2_options // DADA2 
    decontam_options
    taxonomy_options
    database_ch_list
    other_options

    main:

    // Full analysis pipeline (default)
    input_processing_results = WF_INPUT_PROCESSING(
        input_options,
    )

    // Use dynamic metadata if available, otherwise use static path
    metadata_to_use = input_processing_results.metadata_file.ifEmpty(
        metadata_path ? file(metadata_path) : Channel.empty()
    )

    workflow_results = AMPLICON_CORE(
        input_processing_results.read_pairs_ch,
        metadata_to_use,
        input_processing_results.barcode_primer_info_ch,
        primers,
        trimming_options,
        dada2_options,
        decontam_options,
        taxonomy_options,
        database_ch_list,
        other_options,
    )

    emit:
    fastqc_results = workflow_results.fastqc_results
    multiqc_report = workflow_results.multiqc_report
    fastqc_results_after_cutadapt = workflow_results.fastqc_results_after_cutadapt
    multiqc_report_after_cutadapt = workflow_results.multiqc_report_after_cutadapt
    dada2_denoising_stats = workflow_results?.dada2_denoising_stats
    dada2_ft_summary = workflow_results.dada2_ft_summary
    dada2_rep_seqs_summary = workflow_results.dada2_rep_seqs_summary
    dada2_denoising_stats_summary = workflow_results.dada2_denoising_stats_summary
    dada2_base_transition_stats_summary = workflow_results.dada2_base_transition_stats_summary
    decontam_scores_freq = workflow_results.decontam_scores_freq
    decontam_scores_freq_viz = workflow_results.decontam_scores_freq_viz
    decontam_scores_prev = workflow_results.decontam_scores_prev
    decontam_scores_prev_viz = workflow_results.decontam_scores_prev_viz
    taxonomy = workflow_results.taxonomy
    taxonomy_summary = workflow_results.taxonomy_summary
    phylogeny_unrooted = workflow_results.phylogeny_unrooted
    phylogeny_rooted = workflow_results.phylogeny_rooted
    dada2_taxonomy = workflow_results.dada2_taxonomy
    taxonomy_evaluation = workflow_results.taxonomy_evaluation
    taxonomy_comparison = workflow_results.taxonomy_comparison
    gram_staining_evaluation = workflow_results.gram_staining_evaluation
    diversity_results = workflow_results.diversity_results
    feature_table = workflow_results.feature_table
    rep_seqs = workflow_results.rep_seqs
}

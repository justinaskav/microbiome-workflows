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
    database_ch_list
    other_options
    analysis_mode // 'full' or 'qc_only'

    main:

    // Conditional workflow execution based on analysis mode
    if (analysis_mode == 'qc_only') {
        // QC-only analysis
        workflow_results = AMPLICON_QC_ONLY(
            input_options,
            primers,
            trimming_options,
            metadata_path,
        )
    } else {
        // Full analysis pipeline (default)
        input_processing_results = WF_INPUT_PROCESSING(
            input_options,
            primers,
            trimming_options,
        )

        workflow_results = AMPLICON_CORE(
            input_processing_results.read_pairs_ch,
            metadata_path,
            input_processing_results.barcode_primer_info_ch,
            input_processing_results.cutadapt_options,
            dada2_options,
            decontam_options,
            database_ch_list,
            other_options,
        )
    }

    emit:
    // QC results (available in both modes)
    fastqc_results = workflow_results.fastqc_results
    multiqc_report = workflow_results.multiqc_report
    fastqc_results_after_cutadapt = workflow_results.fastqc_results_after_cutadapt
    multiqc_report_after_cutadapt = workflow_results.multiqc_report_after_cutadapt
    
    // Full analysis results (only available in full mode)
    dada2_denoising_stats = workflow_results.dada2_denoising_stats ?: Channel.empty()
    dada2_ft_summary = workflow_results.dada2_ft_summary ?: Channel.empty()
    dada2_rep_seqs_summary = workflow_results.dada2_rep_seqs_summary ?: Channel.empty()
    dada2_denoising_stats_summary = workflow_results.dada2_denoising_stats_summary ?: Channel.empty()
    decontam_scores_freq = workflow_results.decontam_scores_freq ?: Channel.empty()
    decontam_scores_freq_viz = workflow_results.decontam_scores_freq_viz ?: Channel.empty()
    decontam_scores_prev = workflow_results.decontam_scores_prev ?: Channel.empty()
    decontam_scores_prev_viz = workflow_results.decontam_scores_prev_viz ?: Channel.empty()
    taxonomy = workflow_results.taxonomy ?: Channel.empty()
    taxonomy_summary = workflow_results.taxonomy_summary ?: Channel.empty()
    phylogeny_unrooted = workflow_results.phylogeny_unrooted ?: Channel.empty()
    phylogeny_rooted = workflow_results.phylogeny_rooted ?: Channel.empty()
    dada2_taxonomy = workflow_results.dada2_taxonomy ?: Channel.empty()
    taxonomy_evaluation = workflow_results.taxonomy_evaluation ?: Channel.empty()
    taxonomy_comparison = workflow_results.taxonomy_comparison ?: Channel.empty()
    gram_staining_evaluation = workflow_results.gram_staining_evaluation ?: Channel.empty()
    diversity_results = workflow_results.diversity_results ?: Channel.empty()
    feature_table = workflow_results.feature_table ?: Channel.empty()
    rep_seqs = workflow_results.rep_seqs ?: Channel.empty()
}

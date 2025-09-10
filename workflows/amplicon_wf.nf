// Include input processors
include { WF_ZIP_INPUT } from '../modules/input_processors/zip_input_wf.nf'
include { WF_FASTQ_INPUT } from '../modules/input_processors/fastq_input_wf.nf'
include { WF_ACCESSION_INPUT } from '../modules/input_processors/accession_input_wf.nf'
include { WF_QZA_INPUT } from '../modules/input_processors/qza_input_wf.nf'

// Include core processing workflow
include { AMPLICON_CORE } from './amplicon_core_wf.nf'

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

    main:

    // Input processing based on input type
    if (input_options?.zip_filename != null) {
        // ZIP input mode
        input_results = WF_ZIP_INPUT(
            input_options.data_dir,
            input_options.zip_filename,
        )
    }
    else if (input_options?.accession_list) {
        // Accession input mode
        input_results = WF_ACCESSION_INPUT(
            input_options.accession_list,
            input_options.email ?: params.ncbi_email,
            input_options.retries ?: 3,
        )
    }
    else if (input_options?.metadata_qza && input_options?.paired_reads_qza) {
        // QZA input mode
        input_results = WF_QZA_INPUT(
            input_options.metadata_qza,
            input_options.paired_reads_qza,
        )
    }
    else {
        // FASTQ input mode (default)
        input_results = WF_FASTQ_INPUT(
            input_options.data_dir
        )
    }

    // Combine cutadapt options with primer information from batch config
    cutadapt_options = [
        enabled: trimming_options.enabled ?: true,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
        primer_f: primers.forward,
        primer_r: primers.reverse,
        skip_polyg_trim: trimming_options.skip_polyg_trim ?: false,
    ]

    // Run core amplicon processing workflow
    core_results = AMPLICON_CORE(
        input_results.read_pairs_ch,
        metadata_path,
        input_results.barcode_primer_info_ch,
        cutadapt_options,
        dada2_options,
        decontam_options,
        database_ch_list,
        other_options,
    )

    emit:
    fastqc_results = core_results.fastqc_results
    multiqc_report = core_results.multiqc_report
    fastqc_results_after_cutadapt = core_results.fastqc_results_after_cutadapt
    multiqc_report_after_cutadapt = core_results.multiqc_report_after_cutadapt
    dada2_denoising_stats = core_results.dada2_denoising_stats
    dada2_ft_summary = core_results.dada2_ft_summary
    dada2_rep_seqs_summary = core_results.dada2_rep_seqs_summary
    dada2_denoising_stats_summary = core_results.dada2_denoising_stats_summary
    decontam_scores_freq = core_results.decontam_scores_freq
    decontam_scores_freq_viz = core_results.decontam_scores_freq_viz
    decontam_scores_prev = core_results.decontam_scores_prev
    decontam_scores_prev_viz = core_results.decontam_scores_prev_viz
    taxonomy = core_results.taxonomy
    taxonomy_summary = core_results.taxonomy_summary
    phylogeny_unrooted = core_results.phylogeny_unrooted
    phylogeny_rooted = core_results.phylogeny_rooted
    dada2_taxonomy = core_results.dada2_taxonomy
    taxonomy_evaluation = core_results.taxonomy_evaluation
    taxonomy_comparison = core_results.taxonomy_comparison
    diversity_results = core_results.diversity_results
    feature_table = core_results.feature_table
    rep_seqs = core_results.rep_seqs
}

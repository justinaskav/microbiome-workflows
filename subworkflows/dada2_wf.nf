// DADA2 Denoising Subworkflow
// Performs DADA2 denoising on paired-end sequences

nextflow.enable.dsl = 2

include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'

workflow DADA2 {
    take:
    sequences_ch   // Channel of [sample_id, [r1, r2]]
    dada2_options  // Map with trunc_len_f, trunc_len_r, trim_left_f, trim_left_r, n_reads_learn, n_threads, sample_prefix_to_remove

    main:

    dada2_config = [
        trunc_len_f: dada2_options.trunc_len_f,
        trunc_len_r: dada2_options.trunc_len_r,
        trim_left_f: dada2_options.trim_left_f ?: 0,
        trim_left_r: dada2_options.trim_left_r ?: 0,
        max_ee_f: dada2_options.max_ee_f ?: 2,
        max_ee_r: dada2_options.max_ee_r ?: 2,
        n_reads_learn: dada2_options.n_reads_learn ?: params.dada2?.n_reads_learn ?: 1000000,
        n_threads: dada2_options.n_threads ?: params.dada2?.n_threads ?: 8,
        sample_prefix_to_remove: dada2_options.sample_prefix_to_remove ?: params.dada2?.sample_prefix_to_remove,
        binned_quality_scores: dada2_options.binned_quality_scores ?: params.dada2?.binned_quality_scores ?: null,
    ]

    dada2_output_ch = WF_DENOISE_DADA2(
        sequences_ch,
        dada2_config.trunc_len_f,
        dada2_config.trunc_len_r,
        dada2_config.trim_left_f,
        dada2_config.trim_left_r,
        dada2_config.max_ee_f,
        dada2_config.max_ee_r,
        dada2_config.n_reads_learn,
        dada2_config.n_threads,
        dada2_config.sample_prefix_to_remove,
        dada2_config.binned_quality_scores,
    )

    emit:
    feature_table = dada2_output_ch.feature_table
    rep_seqs = dada2_output_ch.representative_sequences
    denoising_stats = dada2_output_ch.denoising_stats
    demux_summary = dada2_output_ch.demux_summary
    feature_table_summary = dada2_output_ch.feature_table_summary
    rep_seqs_summary = dada2_output_ch.rep_seqs_summary
    denoising_stats_summary = dada2_output_ch.denoising_stats_summary
    base_transition_stats = dada2_output_ch.base_transition_stats
}

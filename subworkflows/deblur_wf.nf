// Deblur Denoising Subworkflow
// Performs Deblur denoising on paired-end sequences via read joining

nextflow.enable.dsl = 2

include { WF_DENOISE_DEBLUR } from '../modules/deblur/denoising_deblur_wf.nf'

workflow DEBLUR {
    take:
    sequences_ch    // Channel of [meta, [r1, r2]]
    deblur_options  // Map: trim_length (required), mode, mean_error, indel_prob, min_reads, min_size, jobs_to_start
    reference_seqs  // Channel of path to FeatureData[Sequence] QZA (for mode='other'; Channel.empty() otherwise)

    main:

    deblur_output_ch = WF_DENOISE_DEBLUR(
        sequences_ch,
        deblur_options.trim_length,
        deblur_options.mean_error ?: 0.005,
        deblur_options.indel_prob ?: 0.01,
        deblur_options.min_reads ?: 10,
        deblur_options.min_size ?: 2,
        deblur_options.jobs_to_start ?: params.dada2?.n_threads ?: 1,
        deblur_options.sample_prefix_to_remove ?: params.dada2?.sample_prefix_to_remove,
        deblur_options.mode ?: '16s',
        reference_seqs,
    )

    emit:
    feature_table           = deblur_output_ch.feature_table
    rep_seqs                = deblur_output_ch.representative_sequences
    denoising_stats         = deblur_output_ch.denoising_stats
    demux_summary           = deblur_output_ch.demux_summary
    feature_table_summary   = deblur_output_ch.feature_table_summary
    rep_seqs_summary        = deblur_output_ch.rep_seqs_summary
    denoising_stats_summary = deblur_output_ch.denoising_stats_summary
}

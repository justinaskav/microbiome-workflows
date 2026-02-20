// Trimming Subworkflow
// Runs cutadapt primer/adapter trimming with optional QC after trimming

nextflow.enable.dsl = 2

include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_TRIM } from '../modules/quality_control/quality_control_wf.nf'

workflow TRIMMING {
    take:
    read_pairs_ch          // Channel of [sample_id, [r1, r2]]
    barcode_primer_info_ch // Barcode/primer info (may be empty)
    primers                // Map with forward/reverse primers
    trimming_options       // Map with enabled, use_sampleseq, skip_polyg_trim

    main:

    trimming_config = [
        enabled: trimming_options.enabled ?: true,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
    ]

    if (trimming_config.enabled) {
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
            .set { trimmed_reads_ch }

        // Run QC on trimmed reads
        WF_QUALITY_CONTROL_AFTER_TRIM(trimmed_reads_ch)

        trimmed_sequences = trimming_results.trimmed_sequences
        trimmed_sequences_summary = trimming_results.trimmed_sequences_summary
        fastqc_after_trim = WF_QUALITY_CONTROL_AFTER_TRIM.out.fastqc_results
        multiqc_after_trim = WF_QUALITY_CONTROL_AFTER_TRIM.out.multiqc_report
    } else {
        // No trimming - pass through reads unchanged
        trimmed_reads_ch = read_pairs_ch
        trimmed_sequences = Channel.empty()
        trimmed_sequences_summary = Channel.empty()
        fastqc_after_trim = Channel.empty()
        multiqc_after_trim = Channel.empty()
    }

    emit:
    trimmed_reads = trimmed_reads_ch
    trimmed_sequences = trimmed_sequences
    trimmed_sequences_summary = trimmed_sequences_summary
    fastqc_results_after_trim = fastqc_after_trim
    multiqc_report_after_trim = multiqc_after_trim
}

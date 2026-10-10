// Trimming Subworkflow
// Runs cutadapt primer/adapter trimming with optional QC after trimming

nextflow.enable.dsl = 2

include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_TRIM } from '../modules/quality_control/quality_control_wf.nf'

workflow TRIMMING {
    take:
    read_pairs_ch // Channel of [sample_id, [r1, r2]]
    barcode_primer_info_ch // Barcode/primer info (may be empty)
    trimming_options // Map with enabled, use_sampleseq, nextseq_trim, skip_qc, skip_merged_artifact

    main:

    trimming_config = [
        enabled: trimming_options.enabled,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
        nextseq_trim: trimming_options.nextseq_trim ?: null,
        skip_qc: trimming_options.skip_qc ?: false,
        skip_merged_artifact: trimming_options.skip_merged_artifact ?: false,
    ]

    println "Trimming configuration in trimming_wf.nf: ${trimming_config}"

    if (trimming_config.enabled) {
        // TODO: In use_sampleseq mode, primers come from barcode_primer_info_ch per sample.
        // For non-sampleseq mode, primers would need to be passed explicitly from the batch config.
        def primers = trimming_options.primers ?: [:]
        trimming_results = WF_CUTADAPT(
            read_pairs_ch,
            barcode_primer_info_ch,
            trimming_config.use_sampleseq,
            trimming_config.nextseq_trim,
            primers.forward ?: null,
            primers.reverse ?: null,
            trimming_config.skip_merged_artifact
        )

        // Process exported fastqs for QC after cutadapt
        trimming_results.exported_fastqs
            .map { sample_name, fastq1, fastq2 ->
                return tuple(sample_name, [fastq1, fastq2])
            }
            .set { trimmed_reads_ch }

        trimmed_sequences = trimming_results.trimmed_sequences
        trimmed_sequences_summary = trimming_results.trimmed_sequences_summary
        trimming_stats_ch = trimming_results.trimming_stats

        // Run QC on trimmed reads (optional)
        if (!trimming_config.skip_qc) {
            WF_QUALITY_CONTROL_AFTER_TRIM(trimmed_reads_ch)
            fastqc_after_trim = WF_QUALITY_CONTROL_AFTER_TRIM.out.fastqc_results
            multiqc_after_trim = WF_QUALITY_CONTROL_AFTER_TRIM.out.multiqc_report
        }
        else {
            fastqc_after_trim = Channel.empty()
            multiqc_after_trim = Channel.empty()
        }
    }
    else {
        // No trimming - pass through reads unchanged
        trimmed_reads_ch = read_pairs_ch
        trimmed_sequences = Channel.empty()
        trimmed_sequences_summary = Channel.empty()
        trimming_stats_ch = Channel.empty()
        fastqc_after_trim = Channel.empty()
        multiqc_after_trim = Channel.empty()
    }

    emit:
    trimmed_reads = trimmed_reads_ch
    trimmed_sequences = trimmed_sequences
    trimmed_sequences_summary = trimmed_sequences_summary
    trimming_stats = trimming_stats_ch
    fastqc_results_after_trim = fastqc_after_trim
    multiqc_report_after_trim = multiqc_after_trim
}

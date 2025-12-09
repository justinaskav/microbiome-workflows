// Amplicon Trimming-Only Workflow
// Runs input processing + QC + cutadapt, emits trimmed reads for later merging

nextflow.enable.dsl = 2

include { WF_INPUT_PROCESSING } from './amplicon_input_processing_wf.nf'
include { WF_QUALITY_CONTROL ; WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_CUTADAPT } from '../modules/quality_control/quality_control_wf.nf'
include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'

workflow AMPLICON_TRIMMING_ONLY {
    take:
    name
    metadata_path
    input_options
    primers
    trimming_options

    main:

    input_processing_results = WF_INPUT_PROCESSING(input_options)

    metadata_to_use = input_processing_results.metadata_file.ifEmpty(
        metadata_path ? file(metadata_path) : Channel.empty()
    )

    trimming_config = [
        enabled: trimming_options.enabled ?: true,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
        skip_polyg_trim: trimming_options.skip_polyg_trim ?: false,
    ]

    // QC on original reads
    WF_QUALITY_CONTROL(input_processing_results.read_pairs_ch)

    if (trimming_config.enabled) {
        trimming_results = WF_CUTADAPT(
            input_processing_results.read_pairs_ch,
            input_processing_results.barcode_primer_info_ch,
            trimming_config,
            primers
        )

        trimmed_reads_ch = trimming_results.exported_fastqs
            .map { sample_name, fastq1, fastq2 -> tuple(sample_name, [fastq1, fastq2]) }

        WF_QUALITY_CONTROL_AFTER_CUTADAPT(trimmed_reads_ch)

        fastqc_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.fastqc_results
        multiqc_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.multiqc_report
        trimmed_qza = trimming_results.trimmed_sequences
    } else {
        trimmed_reads_ch = input_processing_results.read_pairs_ch
        fastqc_after_cutadapt = Channel.empty()
        multiqc_after_cutadapt = Channel.empty()
        trimmed_qza = Channel.empty()
    }

    emit:
    trimmed_reads = trimmed_reads_ch
    trimmed_qza = trimmed_qza
    metadata_file = metadata_to_use
    barcode_primer_info = input_processing_results.barcode_primer_info_ch
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    fastqc_results_after_cutadapt = fastqc_after_cutadapt
    multiqc_report_after_cutadapt = multiqc_after_cutadapt
}

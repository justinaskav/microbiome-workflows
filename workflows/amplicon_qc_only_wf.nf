// QC-Only Amplicon Workflow
// Performs FastQC and MultiQC analysis only, with optional cutadapt preprocessing

nextflow.enable.dsl = 2

// Include required workflows and modules
include { WF_QUALITY_CONTROL ; WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_CUTADAPT } from '../modules/quality_control/quality_control_wf.nf'
include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_INPUT_PROCESSING } from './amplicon_input_processing_wf.nf'

// QC-only workflow for FastQC/MultiQC analysis
workflow AMPLICON_QC_ONLY {
    take:
    input_options
    primers
    trimming_options
    metadata_path

    main:

    // Process inputs using the shared input processing workflow
    input_processing_results = WF_INPUT_PROCESSING(
        input_options,
        primers,
        trimming_options,
    )

    read_pairs_ch = input_processing_results.read_pairs_ch
    barcode_primer_info_ch = input_processing_results.barcode_primer_info_ch
    cutadapt_options = input_processing_results.cutadapt_options

    // Run quality control on original reads
    WF_QUALITY_CONTROL(read_pairs_ch)

    // Conditionally run cutadapt if enabled
    if (cutadapt_options.enabled) {
        trimming_results = WF_CUTADAPT(
            read_pairs_ch,
            barcode_primer_info_ch,
            cutadapt_options.enabled,
            cutadapt_options.primer_f,
            cutadapt_options.primer_r,
            cutadapt_options.skip_polyg_trim,
        )

        // Process exported fastqs for QC after cutadapt
        trimming_results.exported_fastqs
            .map { sample_name, fastq1, fastq2 ->
                return tuple(sample_name, [fastq1, fastq2])
            }
            .set { exported_fastqs_pairs_ch }

        // Run QC on the exported fastqs after cutadapt
        WF_QUALITY_CONTROL_AFTER_CUTADAPT(exported_fastqs_pairs_ch)
        
        fastqc_results_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.fastqc_results
        multiqc_report_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.multiqc_report
    } else {
        // No cutadapt - create empty channels
        fastqc_results_after_cutadapt = Channel.empty()
        multiqc_report_after_cutadapt = Channel.empty()
    }

    emit:
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    fastqc_results_after_cutadapt = fastqc_results_after_cutadapt
    multiqc_report_after_cutadapt = multiqc_report_after_cutadapt
}
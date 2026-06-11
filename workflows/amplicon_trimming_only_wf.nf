// Amplicon Trimming-Only Workflow
// Runs input processing + QC + cutadapt, emits trimmed reads for later merging
// Also serves as the de facto QC workflow (provides QC before and after trimming)

nextflow.enable.dsl = 2

// Include subworkflows
include { INPUT_PROCESSING } from '../subworkflows/input_processing_wf.nf'
include { QC } from '../subworkflows/qc_wf.nf'
include { HOST_ALIGNMENT } from '../subworkflows/host_alignment_wf.nf'
include { TRIMMING } from '../subworkflows/trimming_wf.nf'

workflow AMPLICON_TRIMMING_ONLY {
    take:
    name
    metadata_path
    input_options
    trimming_options
    host_options
    host_index_ch

    main:

    // Step 1: Input Processing
    input_processing_results = INPUT_PROCESSING(input_options)

    // Use dynamic metadata if available, otherwise use static path
    metadata_to_use = input_processing_results.metadata_file.ifEmpty(
        metadata_path ? file(metadata_path) : Channel.empty()
    )

    // Step 2: Quality Control on original reads
    qc_results = QC(input_processing_results.read_pairs_ch)

    // Step 3: Trimming
    if (trimming_options?.enabled) {
        trimming_results = TRIMMING(
            input_processing_results.read_pairs_ch,
            input_processing_results.barcode_primer_info_ch,
            trimming_options
        )

        reads_after_trim = trimming_results.trimmed_reads
        trimmed_qza = trimming_results.trimmed_sequences
        fastqc_after_cutadapt = trimming_results.fastqc_results_after_trim
        multiqc_after_cutadapt = trimming_results.multiqc_report_after_trim
    } else {
        reads_after_trim = input_processing_results.read_pairs_ch
        trimmed_qza = Channel.empty()
        fastqc_after_cutadapt = Channel.empty()
        multiqc_after_cutadapt = Channel.empty()
    }

    // Step 3b: Host-alignment QC + optional removal (after cutadapt so bowtie2
    // sees primer-free inserts; default end-to-end mode then aligns correctly).
    host_results = HOST_ALIGNMENT(
        reads_after_trim,
        host_options,
        host_index_ch,
    )
    trimmed_reads_ch = host_results.reads_out

    emit:
    trimmed_reads = trimmed_reads_ch
    trimmed_qza = trimmed_qza
    metadata_file = metadata_to_use
    barcode_primer_info = input_processing_results.barcode_primer_info_ch
    fastqc_results = qc_results.fastqc_results
    multiqc_report = qc_results.multiqc_report
    host_alignment_multiqc = host_results.multiqc_report
    host_alignment_logs = host_results.logs
    fastqc_results_after_cutadapt = fastqc_after_cutadapt
    multiqc_report_after_cutadapt = multiqc_after_cutadapt
}

// Include input processors
include { WF_ZIP_INPUT } from '../modules/input_processors/zip_input_wf.nf'
include { WF_FASTQ_INPUT } from '../modules/input_processors/fastq_input_wf.nf'
include { WF_ACCESSION_INPUT } from '../modules/input_processors/accession_input_wf.nf'
include { WF_QZA_INPUT } from '../modules/input_processors/qza_input_wf.nf'

// Input processing workflow - handles all input types and prepares options
workflow WF_INPUT_PROCESSING {
    take:
    input_options

    main:

    println("Input options: ${input_options}")
    
    // Input processing based on input type
    if (input_options?.zip_filename != null) {
        // ZIP input mode
        input_results = WF_ZIP_INPUT(
            input_options.data_dir,
            input_options.zip_filename,
        )
    }
    else if (input_options?.accession_ids) {
        // Accession input mode
        input_results = WF_ACCESSION_INPUT(
            input_options.accession_ids,
            input_options.email ?: params.ncbi_email,
            input_options?.retries ?: 3,
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
        // Pass full input_options to allow fastq_pattern customization
        input_results = WF_FASTQ_INPUT(input_options)
    }

    emit:
    read_pairs_ch = input_results.read_pairs_ch
    barcode_primer_info_ch = input_results.barcode_primer_info_ch
    metadata_file = input_results.metadata_tsv ?: Channel.empty()
    input_type = input_results.input_type
}

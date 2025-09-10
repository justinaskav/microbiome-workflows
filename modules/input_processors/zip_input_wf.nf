// ZIP Input Processing Workflow
// Processes zipped Novogene data and outputs standardized fastq pairs + metadata

nextflow.enable.dsl = 2

include { WF_UNZIP_DATA } from '../data_processing/unzip_data_wf.nf'
include { WF_PARSE_NOVOGENE_SAMPLESEQ } from '../novogene/sampleseq_wf.nf'

workflow WF_ZIP_INPUT {
    take:
    data_dir                               // Data directory path
    zip_filename                          // Zip filename (without .zip extension)
    // batch_samples_with_additional_barcodes // Samples with additional barcodes for this batch
    // sample_prefix_to_remove               // Optional prefix to remove from sample IDs

    main:
    // Extract data from zip file

    unzip_results = WF_UNZIP_DATA(data_dir, zip_filename)
    fastq_files_ch = unzip_results.fastq_files
    sampleseq_file_ch = unzip_results.sampleseq_file

    // Process fastq files into read pairs
    fastq_files_ch
    .flatMap { list_of_files ->
        list_of_files.collect { file ->
            def file_name = file.getName()
            def sample_name = file_name.replaceAll(/_(?:S[0-9]+_)?L001_R[12]_001\.fastq\.gz$/, "") // Handle optional S number
            return tuple(sample_name, file)
        }
    }
    .groupTuple()
    .set { read_pairs_ch }

    // Parse SampleSeq info for barcode/primer information
    barcode_primer_info_ch = WF_PARSE_NOVOGENE_SAMPLESEQ(
        sampleseq_file_ch,
        // batch_samples_with_additional_barcodes
    ).barcode_primer_info

    emit:
    // Standardized outputs
    read_pairs_ch = read_pairs_ch
    barcode_primer_info_ch = barcode_primer_info_ch
    sampleseq_file_ch = sampleseq_file_ch
    input_type = Channel.value("zip")
    // Cutadapt options for ZIP mode - use SampleSeq info
}

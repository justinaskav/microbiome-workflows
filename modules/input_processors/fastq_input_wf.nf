// FASTQ Input Processing Workflow  
// Processes direct fastq files and outputs standardized fastq pairs + metadata

nextflow.enable.dsl = 2

workflow WF_FASTQ_INPUT {
    take:
    data_dir                    // Data directory path

    main:
    // Data is already in standard format - directly read fastq files
    // Look for standard paired-end fastq files
    read_pairs_ch = Channel.fromFilePairs("${data_dir}/*_L001_{R1,R2}_001.fastq.gz", size: 2)

    // Trim _S part from sample names if present
    read_pairs_ch = read_pairs_ch.map { sample_name, files ->
        def trimmed_sample_name = sample_name.replaceAll(/_S\d+$/, "") // Remove any _S<number> suffix
        return tuple(trimmed_sample_name, files)
    }

    // Create empty sampleseq file channel since no SampleSeq_info.xls is needed
    sampleseq_file_ch = Channel.empty()

    // Create empty barcode_primer_info channel since no SampleSeq parsing needed
    barcode_primer_info_ch = Channel.empty()

    emit:
    // Standardized outputs
    read_pairs_ch = read_pairs_ch
    barcode_primer_info_ch = barcode_primer_info_ch
    sampleseq_file_ch = sampleseq_file_ch
    input_type = Channel.value("fastq")
}

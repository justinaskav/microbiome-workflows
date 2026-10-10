// FASTQ Input Processing Workflow
// Processes direct fastq files and outputs standardized fastq pairs + metadata

nextflow.enable.dsl = 2

workflow WF_FASTQ_INPUT {
    take:
    data_dir                    // Data directory path or map with data_dir and fastq_pattern

    main:
    // Handle both string (legacy) and map input
    def dir_path = data_dir instanceof Map ? data_dir.data_dir : data_dir
    def custom_pattern = data_dir instanceof Map ? data_dir.fastq_pattern : null

    // Use custom pattern, params pattern, or try default patterns
    def pattern = custom_pattern ?: params.fastq_pattern ?: '*_{R1,R2}_001.fastq.gz'

    read_pairs_ch = Channel.fromFilePairs("${dir_path}/${pattern}", size: 2)

    // Trim _S and _L001 parts from sample names if present for cleaner IDs
    read_pairs_ch = read_pairs_ch.map { sample_name, files ->
        def trimmed_sample_name = sample_name
            .replaceAll(/_S\d+$/, "")      // Remove _S<number> suffix
            .replaceAll(/_L\d{3}$/, "")    // Remove _L001 etc suffix
        return tuple(trimmed_sample_name, files)
    }

    // Create empty sampleseq file channel since no SampleSeq_info.xls is needed
    sampleseq_file_ch = Channel.empty()

    // Per-sample barcode/primer table (same columns PARSE_SAMPLESEQ_INFO emits), for
    // barcoded reads whose SampleSeq_info.xls is not shipped with the FASTQs (e.g. ENA).
    def barcode_file = data_dir instanceof Map ? data_dir.barcode_primer_file : null
    barcode_primer_info_ch = barcode_file
        ? Channel.fromPath(barcode_file, checkIfExists: true)
        : Channel.empty()

    emit:
    // Standardized outputs
    read_pairs_ch = read_pairs_ch
    barcode_primer_info_ch = barcode_primer_info_ch
    sampleseq_file_ch = sampleseq_file_ch
    metadata_tsv = Channel.empty()  // No metadata for FASTQ input
    input_type = Channel.value("fastq")
}

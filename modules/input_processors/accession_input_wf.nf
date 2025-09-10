// Accession Input Processing Workflow
// Fetches data via q2-fondue from accession IDs and outputs standardized fastq pairs + metadata

nextflow.enable.dsl = 2

include { WF_QZA_INPUT } from './qza_input_wf.nf'

// Process to fetch data using q2-fondue
process FETCH_ACCESSIONS_FONDUE {
    tag "fetch_accessions_fondue"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"

    input:
    path accession_ids_qza  // QZA file containing accession IDs
    val email              // Email for NCBI/SRA access
    val retries            // Number of retry attempts

    output:
    path "fondue-output/paired-reads.qza", emit: paired_reads_qza
    path "fondue-output/metadata.qza", emit: metadata_qza
    path "fondue-output", emit: output_dir

    script:
    """
    # Create output directory
    mkdir -p fondue-output

    # Run q2-fondue to fetch data
    qiime fondue get-all \\
        --i-accession-ids ${accession_ids_qza} \\
        --p-email ${email} \\
        --p-retries ${retries} \\
        --verbose \\
        --output-dir fondue-output

    # Verify outputs exist
    if [ ! -f "fondue-output/paired-reads.qza" ]; then
        echo "Error: paired-reads.qza not found in fondue output"
        exit 1
    fi
    
    if [ ! -f "fondue-output/metadata.qza" ]; then
        echo "Error: metadata.qza not found in fondue output"
        exit 1
    fi
    """
}

// Process to validate fondue output structure
process VALIDATE_FONDUE_OUTPUT {
    tag "validate_fondue_output"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"

    input:
    path fondue_output_dir
    path paired_reads_qza
    path metadata_qza

    output:
    path paired_reads_qza, emit: validated_paired_reads
    path metadata_qza, emit: validated_metadata
    stdout emit: validation_log

    script:
    """
    echo "Validating fondue output..."
    
    # Check paired-reads.qza type
    qiime tools peek ${paired_reads_qza} | grep -E "Type:|UUID:" || echo "Could not validate paired-reads.qza type"
    
    # Check metadata.qza type  
    qiime tools peek ${metadata_qza} | grep -E "Type:|UUID:" || echo "Could not validate metadata.qza type"
    
    echo "Fondue output validation complete"
    """
}

workflow WF_ACCESSION_INPUT {
    take:
    accession_ids_qza       // QZA file containing accession IDs
    email                   // Email for NCBI/SRA access  
    retries                 // Number of retry attempts (default: 3)

    main:
    // Fetch data using fondue
    fondue_results = FETCH_ACCESSIONS_FONDUE(
        accession_ids_qza,
        email,
        retries ?: 3
    )
    
    // Validate fondue output
    validated_results = VALIDATE_FONDUE_OUTPUT(
        fondue_results.output_dir,
        fondue_results.paired_reads_qza,
        fondue_results.metadata_qza
    )
    
    // Process the downloaded QZA/QZV files using the QZA input workflow
    qza_results = WF_QZA_INPUT(
        validated_results.validated_metadata,
        validated_results.validated_paired_reads
    )

    emit:
    // Standardized outputs matching other input processors
    read_pairs_ch = qza_results.read_pairs_ch
    barcode_primer_info_ch = qza_results.barcode_primer_info_ch
    sampleseq_file_ch = qza_results.sampleseq_file_ch
    input_type = Channel.value("accession")
}

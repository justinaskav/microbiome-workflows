// QZA Input Processing Workflow
// Processes QIIME2 QZA artifacts and outputs standardized fastq pairs + metadata

nextflow.enable.dsl = 2

// Process to export QZA metadata to TSV format
process EXPORT_QZA_METADATA {
    tag "export_qza_metadata"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"

    input:
    path metadata_qza

    output:
    path "metadata.tsv", emit: metadata_tsv
    stdout emit: export_log

    script:
    """
    echo "Exporting QZA metadata to TSV format..."
    
    # Export metadata QZA to TSV
    qiime tools export \\
        --input-path ${metadata_qza} \\
        --output-path exported_metadata
    
    # Find the metadata file (usually named metadata.tsv)
    if [ -f "exported_metadata/metadata.tsv" ]; then
        cp exported_metadata/metadata.tsv metadata.tsv
    elif [ -f "exported_metadata/sample-metadata.tsv" ]; then
        cp exported_metadata/sample-metadata.tsv metadata.tsv
    else
        # Find any TSV file in the exported directory
        tsv_file=\$(find exported_metadata -name "*.tsv" | head -1)
        if [ -n "\$tsv_file" ]; then
            cp "\$tsv_file" metadata.tsv
        else
            echo "Error: No TSV metadata file found in exported QZA"
            exit 1
        fi
    fi
    
    echo "QZA metadata export complete"
    """
}

// Process to export QZA paired reads to FASTQ format
process EXPORT_QZA_PAIRED_READS {
    tag "export_qza_paired_reads"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"

    input:
    path paired_reads_qza

    output:
    path "exported_fastq/*.fastq.gz", emit: fastq_files
    stdout emit: export_log

    script:
    """
    echo "Exporting QZA paired reads to FASTQ format..."
    
    # Export paired reads QZA to FASTQ
    qiime tools export \\
        --input-path ${paired_reads_qza} \\
        --output-path exported_fastq
    
    # Verify FASTQ files were created
    fastq_count=\$(find exported_fastq -name "*.fastq.gz" | wc -l)
    if [ \$fastq_count -eq 0 ]; then
        echo "Error: No FASTQ files found in exported QZA"
        exit 1
    fi
    
    echo "Found \$fastq_count FASTQ files"
    echo "QZA paired reads export complete"
    """
}

// Process to organize FASTQ files into paired-end format
process ORGANIZE_FASTQ_PAIRS {
    tag "organize_fastq_pairs"

    input:
    path "fastq_files/*"

    output:
    tuple val("organized"), path("organized_fastq/*.fastq.gz"), emit: organized_files

    script:
    """
    mkdir -p organized_fastq
    
    # Copy all FASTQ files to organized directory
    cp fastq_files/*.fastq.gz organized_fastq/
    
    # List files for debugging
    echo "Organized FASTQ files:"
    ls -la organized_fastq/
    """
}

workflow WF_QZA_INPUT {
    take:
    metadata_qza           // QZA file containing sample metadata
    paired_reads_qza       // QZA file containing paired-end reads

    main:
    // Export QZA metadata to TSV format
    metadata_results = EXPORT_QZA_METADATA(metadata_qza)
    
    // Export QZA paired reads to FASTQ format
    fastq_results = EXPORT_QZA_PAIRED_READS(paired_reads_qza)
    
    // Organize FASTQ files
    organized_results = ORGANIZE_FASTQ_PAIRS(fastq_results.fastq_files.collect())
    
    // Process FASTQ files into read pairs
    // Extract sample names from FASTQ filenames and group by sample
    read_pairs_ch = organized_results.organized_files
        .flatMap { tuple ->
            def organized_files = tuple[1]
            return organized_files.collect { file ->
                def file_name = file.getName()
                // Extract sample name from various FASTQ naming patterns
                def sample_name = file_name.replaceAll(/(_L001)?_R[12](_001)?\\.fastq\\.gz$/, "")
                    .replaceAll(/_[0-9]+_[0-9]+_forward|_[0-9]+_[0-9]+_reverse/, "")
                    .replaceAll(/_forward|_reverse/, "")
                return tuple(sample_name, file)
            }
        }
        .groupTuple()
        .map { sample_name, files ->
            // Sort files to ensure R1 comes before R2
            def sorted_files = files.sort { it.getName() }
            return tuple(sample_name, sorted_files)
        }

    // Create empty sampleseq file channel since no SampleSeq parsing needed for QZA input
    sampleseq_file_ch = Channel.empty()

    // Create empty barcode_primer_info channel since QZA input doesn't contain barcode info
    barcode_primer_info_ch = Channel.empty()

    emit:
    // Standardized outputs matching other input processors
    read_pairs_ch = read_pairs_ch
    barcode_primer_info_ch = barcode_primer_info_ch
    sampleseq_file_ch = sampleseq_file_ch
    input_type = Channel.value("qza")
    metadata_tsv = metadata_results.metadata_tsv
}
// Process already unzipped Novogene data
// For cases where data is already extracted and doesn't need unzipping

process PROCESS_UNZIPPED_NOVOGENE_DATA {
    tag "Process unzipped data from ${data_dir}"

    input:
    val data_dir
    val batch_name

    output:
    path "renamed/*", emit: fastq_files
    path "SampleSeq_info.xls", emit: sampleseq_file, optional: true

    script:
    """
    # Create directory for renamed files
    mkdir -p renamed

    # Look for SampleSeq_info.xls in the data directory
    if [ -f "${data_dir}/SampleSeq_info.xls" ]; then
        cp "${data_dir}/SampleSeq_info.xls" .
    else
        echo "No SampleSeq_info.xls found in ${data_dir}"
    fi

    # Find all fastq.gz files that match the pattern *.raw_*.fastq.gz
    for file in ${data_dir}/*.raw_*.fastq.gz; do
        if [ -f "\$file" ]; then
            # Extract the sample name (everything before .raw_)
            sample=\$(basename "\$file" | sed 's/\\.raw_.*//g')
            
            # Determine if it's R1 or R2 from the filename
            if [[ \$file == *"raw_1"* ]]; then
                cp "\$file" "renamed/\${sample}_L001_R1_001.fastq.gz"
            elif [[ \$file == *"raw_2"* ]]; then
                cp "\$file" "renamed/\${sample}_L001_R2_001.fastq.gz"
            fi
        fi
    done

    # Check if we found any files
    if [ ! "\$(ls -A renamed)" ]; then
        echo "ERROR: No raw fastq files found in ${data_dir}"
        echo "Expected files matching pattern: *.raw_*.fastq.gz"
        exit 1
    fi
    """
}

// Workflow for processing already unzipped data
workflow WF_PROCESS_UNZIPPED_DATA {
    take:
    data_dir
    batch_name

    main:
    process_results = PROCESS_UNZIPPED_NOVOGENE_DATA(data_dir, batch_name)

    emit:
    fastq_files = process_results.fastq_files
    sampleseq_file = process_results.sampleseq_file
}

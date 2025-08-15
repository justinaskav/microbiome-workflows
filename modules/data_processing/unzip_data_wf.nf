// Extract data from a single zip file

process EXTRACT_NOVOGENE_ZIP_DATA {
    tag "Extract data from ${zip_filename}"

    input:
    path path_to_zip
    val zip_filename

    output:
    path "renamed/*", emit: fastq_files
    path "SampleSeq_info.xls", emit: sampleseq_file, optional: true

    script:
    """
    # First unzip the files
    mkdir -p unzipped
    unzip -j -o ${path_to_zip} '*.raw*.gz' -d unzipped
    
    # Also extract SampleSeq_info.xls if it exists
    unzip -j -o ${path_to_zip} '*SampleSeq_info.xls' -d . || echo "No SampleSeq_info.xls found"

    # Create directory for renamed files
    mkdir -p renamed

    # Loop through the files and rename them
    for file in unzipped/*.raw_*.fastq.gz; do
        # Extract the sample name (everything before .raw_)
        sample=\$(basename "\$file" | sed 's/\\.raw_.*//g')
        
        # Determine if it's R1 or R2 from the filename
        if [[ \$file == *"raw_1"* ]]; then
            mv "\$file" "renamed/\${sample}_L001_R1_001.fastq.gz"
        elif [[ \$file == *"raw_2"* ]]; then
            mv "\$file" "renamed/\${sample}_L001_R2_001.fastq.gz"
        fi
    done

    # Clean up
    rm -rf unzipped
    """
}

// Simple workflow to extract data from a single zip file
workflow WF_UNZIP_DATA {
    take:
    path_to_zip
    zip_filename

    main:
    extract_results = EXTRACT_NOVOGENE_ZIP_DATA(path_to_zip, zip_filename)

    emit:
    fastq_files = extract_results.fastq_files
    sampleseq_file = extract_results.sampleseq_file
}

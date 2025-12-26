// Merge trimmed reads from multiple batches before DADA2
// This allows DADA2 to learn error models from combined data

nextflow.enable.dsl = 2

// Simple process to collect and pass through trimmed reads
// The actual "merge" is just combining channels - DADA2 handles them together
process COLLECT_TRIMMED_READS {
    tag "Collect trimmed reads"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.fastq.gz"), emit: collected_reads

    script:
    """
    # Rename to consistent format to avoid conflicts
    cp -L ${reads[0]} ${sample_id}_R1.fastq.gz
    cp -L ${reads[1]} ${sample_id}_R2.fastq.gz
    """
}

// Merge metadata files from multiple batches
process MERGE_METADATA {
    tag "Merge metadata"

    input:
    path metadata_files, stageAs: "metadata_*.tsv"

    output:
    path "merged_metadata.tsv", emit: merged_metadata

    script:
    """
    # QIIME2 metadata has 2 header lines: column names and #q2:types directive
    # Take both header lines from first file, then data from all files (skipping both headers)
    first_file=\$(ls metadata_*.tsv | head -1)
    head -2 "\$first_file" > merged_metadata.tsv
    for f in metadata_*.tsv; do
        tail -n +3 "\$f" >> merged_metadata.tsv
    done
    """
}

workflow WF_MERGE_TRIMMED_READS {
    take:
    trimmed_reads_ch  // Mixed channel from multiple batches: [sample_id, [r1, r2]]

    main:
    // Collect reads with consistent naming
    collected_ch = COLLECT_TRIMMED_READS(trimmed_reads_ch)

    emit:
    merged_reads = collected_ch.collected_reads
}

workflow WF_MERGE_BATCH_METADATA {
    take:
    metadata_files_ch  // Channel of metadata file paths

    main:
    MERGE_METADATA(metadata_files_ch.collect())

    emit:
    merged_metadata = MERGE_METADATA.out.merged_metadata
}

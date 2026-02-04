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
// Handles QIIME2 format with #q2:types row and aligns columns by name
process MERGE_METADATA {
    tag "Merge metadata"
    conda "conda-forge::pandas"

    input:
    path metadata_files, stageAs: "metadata_*.tsv"

    output:
    path "merged_metadata.tsv", emit: merged_metadata

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from pathlib import Path

    files = sorted(Path('.').glob('metadata_*.tsv'))
    dfs = []

    for f in files:
        # Read TSV, skip comment lines starting with #
        df = pd.read_csv(f, sep='\\t', comment='#', dtype=str)
        dfs.append(df)

    # Merge all dataframes, aligning by column name
    merged = pd.concat(dfs, ignore_index=True)

    # Fill missing DecontaminationStatus with 'No' (for batches without controls)
    if 'DecontaminationStatus' in merged.columns:
        merged['DecontaminationStatus'] = merged['DecontaminationStatus'].fillna('No')

    # Write output with QIIME2 types row
    cols = merged.columns.tolist()
    types = ['categorical'] * len(cols)

    with open('merged_metadata.tsv', 'w') as out:
        out.write('\\t'.join(cols) + '\\n')
        out.write('#q2:types\\t' + '\\t'.join(types[1:]) + '\\n')
        merged.to_csv(out, sep='\\t', index=False, header=False)
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

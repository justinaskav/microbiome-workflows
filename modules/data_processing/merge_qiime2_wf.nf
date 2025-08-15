// QIIME2 merging processes for multiple batch processing

// Process to merge feature tables from multiple batches
process MERGE_FEATURE_TABLES {
    tag "Merge feature tables"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"
    publishDir "${params.outputDir}/merged", mode: 'copy'

    input:
    path tables, stageAs: "table_*.qza"

    output:
    path "final_table.qza", emit: merged_table

    script:
    """
    # Merge tables
    qiime feature-table merge \\
      ${tables.collect { file -> "--i-tables ${file}" }.join(' ')} \\
      --o-merged-table final_table.qza
    """
}

// Process to merge representative sequences
process MERGE_SEQUENCES {
    tag "Merge sequences"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"
    publishDir "${params.outputDir}/merged", mode: 'copy'

    input:
    path sequences, stageAs: "rep-seqs_*.qza"

    output:
    path "final_seqs.qza", emit: merged_seqs

    script:
    """
    # Merge sequences
    qiime feature-table merge-seqs \\
      ${sequences.collect { file -> "--i-data ${file}" }.join(' ')} \\
      --o-merged-data final_seqs.qza
    """
}

// Process to merge taxonomies
process MERGE_TAXONOMIES {
    tag "Merge taxonomies"
    conda "${params.qiime2_conda_env ?: 'qiime2-amplicon'}"
    publishDir "${params.outputDir}/merged", mode: 'copy'

    input:
    path taxonomies, stageAs: "taxonomy_*.qza"

    output:
    path "final_taxonomy.qza", emit: merged_taxonomy

    script:
    """
    # Merge taxonomies
    qiime feature-table merge-taxa \\
      ${taxonomies.collect { file -> "--i-data ${file}" }.join(' ')} \\
      --o-merged-data final_taxonomy.qza
    """
}
// Taxonomy visualization workflows for QIIME2
// This module contains reusable processes for visualizing taxonomy data

// Process to visualize taxonomic assignments as a barplot
process VISUALIZE_TAXONOMY_BARPLOT {
    tag "visualize_taxonomy_${taxonomy_qza.getBaseName()}"
    label 'process_low'
    
    conda "${params.qiime2_conda_env}"
    // container (params.qiime2_container ?: "quay.io/qiime2/amplicon:2025.4")

    input:
    path taxonomy_qza        // Taxonomic classification (QZA)
    path feature_table_qza   // Feature table (QZA) for barplot context
    val output_name          // Name for the output file (without extension)
    path metadata            // Metadata file for barplot visualization

    output:
    path("${output_name}.qzv"), emit: barplot_qzv

    script:
    metadata_str = metadata ? "--m-metadata-file ${metadata}" : ""
    """
    qiime taxa barplot \\
        --i-table "${feature_table_qza}" \\
        --i-taxonomy "${taxonomy_qza}" \\
        --o-visualization ${output_name}.qzv \\
        ${metadata_str} \\
        --verbose
    """
}

// A workflow that can be included in other workflows
workflow VISUALIZE_TAXONOMY {
    take:
    taxonomy_qza        // Path to taxonomy artifact (.qza)
    feature_table_qza   // Path to feature table (.qza)
    output_name         // Name for the output visualization
    metadata
    
    main:
    barplot_ch = VISUALIZE_TAXONOMY_BARPLOT(
        taxonomy_qza,
        feature_table_qza,
        output_name,
        metadata
    )
    
    emit:
    barplot = barplot_ch.barplot_qzv
}

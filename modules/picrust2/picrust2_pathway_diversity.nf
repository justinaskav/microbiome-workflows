// PICRUSt2 Pathway Diversity Analysis
// Runs core-metrics (non-phylogenetic) on pathway abundance data

nextflow.enable.dsl = 2

process PICRUST2_PATHWAY_DIVERSITY {
    tag "picrust2_pathway_diversity"
    conda "${params.qiime2_conda_env}"

    input:
    path pathway_abundance_qza
    path metadata_file
    val sampling_depth

    output:
    path "pathway_core_metrics/rarefied_table.qza", emit: rarefied_table
    path "pathway_core_metrics/observed_features_vector.qza", emit: observed_features_vector
    path "pathway_core_metrics/shannon_vector.qza", emit: shannon_vector
    path "pathway_core_metrics/evenness_vector.qza", emit: evenness_vector
    path "pathway_core_metrics/jaccard_distance_matrix.qza", emit: jaccard_distance_matrix
    path "pathway_core_metrics/bray_curtis_distance_matrix.qza", emit: bray_curtis_distance_matrix
    path "pathway_core_metrics/jaccard_pcoa_results.qza", emit: jaccard_pcoa_results
    path "pathway_core_metrics/bray_curtis_pcoa_results.qza", emit: bray_curtis_pcoa_results
    path "pathway_core_metrics/jaccard_emperor.qzv", emit: jaccard_emperor
    path "pathway_core_metrics/bray_curtis_emperor.qzv", emit: bray_curtis_emperor

    script:
    """
    qiime diversity core-metrics \\
        --i-table ${pathway_abundance_qza} \\
        --p-sampling-depth ${sampling_depth} \\
        --m-metadata-file ${metadata_file} \\
        --output-dir pathway_core_metrics \\
        --verbose
    """
}

// Workflow wrapper for pathway diversity
workflow WF_PICRUST2_PATHWAY_DIVERSITY {
    take:
    pathway_abundance_ch    // path to pathway abundance QZA
    metadata_ch             // path to metadata file
    sampling_depth          // integer sampling depth

    main:
    diversity_results = PICRUST2_PATHWAY_DIVERSITY(
        pathway_abundance_ch,
        metadata_ch,
        sampling_depth
    )

    emit:
    rarefied_table = diversity_results.rarefied_table
    observed_features_vector = diversity_results.observed_features_vector
    shannon_vector = diversity_results.shannon_vector
    evenness_vector = diversity_results.evenness_vector
    jaccard_distance_matrix = diversity_results.jaccard_distance_matrix
    bray_curtis_distance_matrix = diversity_results.bray_curtis_distance_matrix
    jaccard_pcoa_results = diversity_results.jaccard_pcoa_results
    bray_curtis_pcoa_results = diversity_results.bray_curtis_pcoa_results
    jaccard_emperor = diversity_results.jaccard_emperor
    bray_curtis_emperor = diversity_results.bray_curtis_emperor
}

// QIIME2 diversity core-metrics-phylogenetic workflow
// This module runs the core diversity metrics phylogenetic analysis

process DIVERSITY_CORE_METRICS_PHYLOGENETIC {
    tag "Core metrics on ${meta.id ?: 'data'}"
    conda params.qiime2_conda_env

    input:
    path feature_table
    path rooted_tree
    path metadata_file
    val sampling_depth
    val meta

    output:
    path "core_metrics_results/rarefied_table.qza", emit: rarefied_table
    path "core_metrics_results/faith_pd_vector.qza", emit: faith_pd_vector
    path "core_metrics_results/observed_features_vector.qza", emit: observed_features_vector
    path "core_metrics_results/shannon_vector.qza", emit: shannon_vector
    path "core_metrics_results/evenness_vector.qza", emit: evenness_vector
    path "core_metrics_results/unweighted_unifrac_distance_matrix.qza", emit: unweighted_unifrac_distance_matrix
    path "core_metrics_results/weighted_unifrac_distance_matrix.qza", emit: weighted_unifrac_distance_matrix
    path "core_metrics_results/jaccard_distance_matrix.qza", emit: jaccard_distance_matrix
    path "core_metrics_results/bray_curtis_distance_matrix.qza", emit: bray_curtis_distance_matrix
    path "core_metrics_results/unweighted_unifrac_pcoa_results.qza", emit: unweighted_unifrac_pcoa_results
    path "core_metrics_results/weighted_unifrac_pcoa_results.qza", emit: weighted_unifrac_pcoa_results
    path "core_metrics_results/jaccard_pcoa_results.qza", emit: jaccard_pcoa_results
    path "core_metrics_results/bray_curtis_pcoa_results.qza", emit: bray_curtis_pcoa_results
    path "core_metrics_results/unweighted_unifrac_emperor.qzv", emit: unweighted_unifrac_emperor
    path "core_metrics_results/weighted_unifrac_emperor.qzv", emit: weighted_unifrac_emperor
    path "core_metrics_results/jaccard_emperor.qzv", emit: jaccard_emperor
    path "core_metrics_results/bray_curtis_emperor.qzv", emit: bray_curtis_emperor

    script:
    """
    qiime diversity core-metrics-phylogenetic \\
        --i-phylogeny ${rooted_tree} \\
        --i-table ${feature_table} \\
        --p-sampling-depth ${sampling_depth} \\
        --m-metadata-file ${metadata_file} \\
        --output-dir core_metrics_results
    """
}

// QIIME2 diversity core-metrics (non-phylogenetic) workflow
// This module runs the core diversity metrics analysis without phylogeny

process DIVERSITY_CORE_METRICS {
    tag "Core metrics on ${meta.id ?: 'data'}"
    conda params.qiime2_conda_env

    input:
    path feature_table
    path metadata_file
    val sampling_depth
    val meta

    output:
    path "core_metrics_results/rarefied_table.qza", emit: rarefied_table
    path "core_metrics_results/observed_features_vector.qza", emit: observed_features_vector
    path "core_metrics_results/shannon_vector.qza", emit: shannon_vector
    path "core_metrics_results/evenness_vector.qza", emit: evenness_vector
    path "core_metrics_results/jaccard_distance_matrix.qza", emit: jaccard_distance_matrix
    path "core_metrics_results/bray_curtis_distance_matrix.qza", emit: bray_curtis_distance_matrix
    path "core_metrics_results/jaccard_pcoa_results.qza", emit: jaccard_pcoa_results
    path "core_metrics_results/bray_curtis_pcoa_results.qza", emit: bray_curtis_pcoa_results
    path "core_metrics_results/jaccard_emperor.qzv", emit: jaccard_emperor
    path "core_metrics_results/bray_curtis_emperor.qzv", emit: bray_curtis_emperor

    script:
    """
    qiime diversity core-metrics \\
        --i-table ${feature_table} \\
        --p-sampling-depth ${sampling_depth} \\
        --m-metadata-file ${metadata_file} \\
        --output-dir core_metrics_results
    """
}

// Workflow to run core metrics and significance tests
workflow WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC {
    take:
    feature_table_ch     // tuple val(meta), path(feature_table)
    rooted_tree_ch       // path rooted_tree
    metadata_ch          // path metadata_file
    sampling_depth       // val sampling_depth
    meta                // val meta

    main:
    // Run core metrics phylogenetic
    core_metrics_results = DIVERSITY_CORE_METRICS_PHYLOGENETIC(
        feature_table_ch,
        rooted_tree_ch,
        metadata_ch,
        sampling_depth,
        meta
    )

    emit:
    // Core metrics results
    rarefied_table = core_metrics_results.rarefied_table
    faith_pd_vector = core_metrics_results.faith_pd_vector
    observed_features_vector = core_metrics_results.observed_features_vector
    shannon_vector = core_metrics_results.shannon_vector
    evenness_vector = core_metrics_results.evenness_vector
    unweighted_unifrac_distance_matrix = core_metrics_results.unweighted_unifrac_distance_matrix
    weighted_unifrac_distance_matrix = core_metrics_results.weighted_unifrac_distance_matrix
    jaccard_distance_matrix = core_metrics_results.jaccard_distance_matrix
    bray_curtis_distance_matrix = core_metrics_results.bray_curtis_distance_matrix
    unweighted_unifrac_pcoa_results = core_metrics_results.unweighted_unifrac_pcoa_results
    weighted_unifrac_pcoa_results = core_metrics_results.weighted_unifrac_pcoa_results
    jaccard_pcoa_results = core_metrics_results.jaccard_pcoa_results
    bray_curtis_pcoa_results = core_metrics_results.bray_curtis_pcoa_results
    unweighted_unifrac_emperor = core_metrics_results.unweighted_unifrac_emperor
    weighted_unifrac_emperor = core_metrics_results.weighted_unifrac_emperor
    jaccard_emperor = core_metrics_results.jaccard_emperor
    bray_curtis_emperor = core_metrics_results.bray_curtis_emperor
}

// Workflow to run non-phylogenetic core metrics (for ITS data)
workflow WF_DIVERSITY_CORE_METRICS {
    take:
    feature_table_ch     // tuple val(meta), path(feature_table)
    metadata_ch          // path metadata_file
    sampling_depth       // val sampling_depth
    meta                // val meta

    main:
    // Run core metrics (non-phylogenetic)
    core_metrics_results = DIVERSITY_CORE_METRICS(
        feature_table_ch,
        metadata_ch,
        sampling_depth,
        meta
    )

    emit:
    // Core metrics results (non-phylogenetic)
    rarefied_table = core_metrics_results.rarefied_table
    observed_features_vector = core_metrics_results.observed_features_vector
    shannon_vector = core_metrics_results.shannon_vector
    evenness_vector = core_metrics_results.evenness_vector
    jaccard_distance_matrix = core_metrics_results.jaccard_distance_matrix
    bray_curtis_distance_matrix = core_metrics_results.bray_curtis_distance_matrix
    jaccard_pcoa_results = core_metrics_results.jaccard_pcoa_results
    bray_curtis_pcoa_results = core_metrics_results.bray_curtis_pcoa_results
    jaccard_emperor = core_metrics_results.jaccard_emperor
    bray_curtis_emperor = core_metrics_results.bray_curtis_emperor
}

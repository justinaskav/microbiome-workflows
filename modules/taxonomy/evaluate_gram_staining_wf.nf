process EVALUATE_GRAM_STAINING {
    conda "${params.qiime2_conda_env}" // Uses the QIIME2 environment
    cache false

    input:
    path taxonomy
    path feature_table
    path metadata_file

    output:
    path "gram_staining_eval_results.qzv", emit: eval_results

    script:
    """
    qiime jk visualize-gram-staining \\
        --i-taxonomy ${taxonomy} \\
        --i-feature-table ${feature_table} \\
        --m-metadata-file ${metadata_file} \\
        --o-visualization gram_staining_eval_results.qzv \\
        --verbose
    """
}

// Workflow for comparing taxonomy artifacts
workflow WF_EVALUATE_GRAM_STAINING {
    take:
    taxonomy         // Channel with taxonomy artifact
    feature_table    // Channel with feature table artifact
    metadata_file    // Channel with metadata file
    
    main:
    // Execute comparison
    eval_ch = EVALUATE_GRAM_STAINING(
        taxonomy,
        feature_table,
        metadata_file
    )
    
    emit:
    gram_staining_eval_results = eval_ch.eval_results
}

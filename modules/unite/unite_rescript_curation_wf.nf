// UNITE database curation workflow using RESCRIPt
// Downloads and curates UNITE database for fungal taxonomic classification
// Based on: https://forum.qiime2.org/t/how-to-train-a-unite-classifier-using-rescript/28285

nextflow.enable.dsl = 2

process UNITE_DOWNLOAD_RESCRIPT {
    tag "Download UNITE database v${unite_version}"
    conda params.qiime2_conda_env

    input:
    val unite_version
    val taxon_group
    val cluster_id

    output:
    path "uniteDB/sequences.qza", emit: sequences
    path "uniteDB/taxonomy.qza", emit: taxonomy
    path "uniteDB/", emit: unite_dir

    script:
    """
    qiime rescript get-unite-data \\
        --p-version ${unite_version} \\
        --p-taxon-group ${taxon_group} \\
        --p-cluster-id ${cluster_id} \\
        --p-no-singletons \\
        --verbose \\
        --output-dir uniteDB
    """
}

process UNITE_FILTER_SEQUENCES {
    tag "Filter UNITE sequences"
    conda params.qiime2_conda_env

    input:
    path sequences
    path taxonomy
    val exclude_terms

    output:
    path "sequences-filtered.qza", emit: filtered_sequences

    script:
    def exclude_list = exclude_terms.join(',')
    """
    qiime taxa filter-seqs \\
        --p-exclude ${exclude_list} \\
        --i-taxonomy ${taxonomy} \\
        --i-sequences ${sequences} \\
        --o-filtered-sequences sequences-filtered.qza
    """
}


process UNITE_EDIT_TAXONOMY {
    tag "Remove SH identifiers from taxonomy"
    conda params.qiime2_conda_env

    input:
    path taxonomy

    output:
    path "taxonomy-no-SH.qza", emit: edited_taxonomy

    script:
    """
    qiime rescript edit-taxonomy \\
        --i-taxonomy ${taxonomy} \\
        --o-edited-taxonomy taxonomy-no-SH.qza \\
        --p-search-strings ';sh__.*' \\
        --p-replacement-strings '' \\
        --p-use-regex
    """
}

process UNITE_TRAIN_CLASSIFIER {
    tag "Train UNITE Naive Bayes classifier"
    conda params.qiime2_conda_env

    input:
    path sequences
    path taxonomy

    output:
    path "unite-classifier.qza", emit: classifier

    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \\
        --i-reference-reads ${sequences} \\
        --i-reference-taxonomy ${taxonomy} \\
        --o-classifier unite-classifier.qza
    """
}

process UNITE_EVALUATE_CLASSIFIER {
    tag "Evaluate UNITE classifier"
    conda params.qiime2_conda_env

    input:
    path sequences
    path taxonomy
    val n_jobs

    output:
    path "classifier-evaluation.qzv", emit: evaluation
    path "predicted-taxonomy.qza", emit: predicted_taxonomy

    script:
    """
    qiime rescript evaluate-fit-classifier \\
        --i-sequences ${sequences} \\
        --i-taxonomy ${taxonomy} \\
        --p-n-jobs ${n_jobs} \\
        --o-classifier temp-classifier.qza \\
        --o-evaluation classifier-evaluation.qzv \\
        --o-observed-taxonomy predicted-taxonomy.qza
    """
}

// Main workflow for UNITE database curation using RESCRIPt
workflow WF_UNITE_RESCRIPT_CURATION {
    take:
    unite_version      // val: UNITE database version
    taxon_group       // val: taxonomic group 
    cluster_id        // val: clustering method
    exclude_terms     // list: terms to exclude during filtering
    n_jobs           // val: number of parallel jobs for evaluation
    enable_evaluation // val: whether to run classifier evaluation

    main:
    // Download UNITE database
    unite_download = UNITE_DOWNLOAD_RESCRIPT(
        unite_version,
        taxon_group,
        cluster_id
    )

    // Remove species hypothesis (SH) identifiers first
    edited_taxonomy = UNITE_EDIT_TAXONOMY(
        unite_download.taxonomy
    )

    // Skip sequence filtering for now - use all sequences with edited taxonomy
    // Filter sequences to remove generic annotations (using edited taxonomy)
    // filtered_seqs = UNITE_FILTER_SEQUENCES(
    //     unite_download.sequences,
    //     edited_taxonomy.edited_taxonomy,
    //     exclude_terms
    // )

    // Train Naive Bayes classifier using all sequences with edited taxonomy
    classifier = UNITE_TRAIN_CLASSIFIER(
        unite_download.sequences,
        edited_taxonomy.edited_taxonomy
    )

    // Optional: Evaluate classifier performance
    if (enable_evaluation) {
        evaluation = UNITE_EVALUATE_CLASSIFIER(
            unite_download.sequences,
            edited_taxonomy.edited_taxonomy,
            n_jobs
        )
        
        classifier_evaluation = evaluation.evaluation
        predicted_taxonomy = evaluation.predicted_taxonomy
    } else {
        classifier_evaluation = Channel.empty()
        predicted_taxonomy = Channel.empty()
    }

    emit:
    // Main outputs
    classifier = classifier.classifier
    sequences = unite_download.sequences
    taxonomy = edited_taxonomy.edited_taxonomy
    
    // Quality control outputs
    evaluation = classifier_evaluation
    predicted_taxonomy = predicted_taxonomy
    
    // Raw database files
    unite_dir = unite_download.unite_dir
}
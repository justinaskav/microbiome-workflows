// Unified classifier creation workflow
// Creates naive bayes classifiers for both SILVA and UNITE databases

process CREATE_NB_CLASSIFIER {
    tag "$meta.id - ${db_type.toUpperCase()} Classifier"
    label 'process_high' 

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(reference_sequences_qza), path(reference_taxonomy_qza)
    val db_type  // 'silva' or 'unite' for naming

    output:
    path("${meta.id}-${db_type}-classifier.qza"), emit: classifier

    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "${reference_sequences_qza}" \
        --i-reference-taxonomy "${reference_taxonomy_qza}" \
        --p-feat-ext--ngram-range '[8, 8]' \
        --o-classifier "${meta.id}-${db_type}-classifier.qza" \
        --verbose
    """
}

process EVALUATE_CLASSIFIER {
    tag "$meta.id - Evaluate ${db_type.toUpperCase()}"
    label 'process_high'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val db_type
    val n_jobs

    output:
    path("${meta.id}-${db_type}-evaluation.qzv"), emit: evaluation
    path("${meta.id}-${db_type}-predicted.qza"), emit: predicted_taxonomy

    when:
    task.ext.run_evaluation == true

    script:
    """
    qiime rescript evaluate-fit-classifier \
        --i-sequences "${sequences_qza}" \
        --i-taxonomy "${taxonomy_qza}" \
        --p-n-jobs ${n_jobs ?: 1} \
        --o-classifier temp-classifier.qza \
        --o-evaluation "${meta.id}-${db_type}-evaluation.qzv" \
        --o-observed-taxonomy "${meta.id}-${db_type}-predicted.qza"
    """
}

include { WF_PROCESS_REF_DB } from './process_reference_db_wf.nf'

workflow CREATE_REF_DB_CLASSIFIER {
    take:
    meta
    db_type              // 'silva' or 'unite'
    db_config            // Database download configuration
    processing_config    // Processing configuration
    classifier_config    // Classifier configuration (evaluation settings, etc.)

    main:
    // Process the reference database
    processed_result = WF_PROCESS_REF_DB(
        meta,
        db_type,
        db_config,
        processing_config
    )
    processed_data_ch = processed_result.processed_data

    // Create the classifier
    classifier_ch = CREATE_NB_CLASSIFIER(
        processed_data_ch,
        db_type
    ).classifier

    // Optional evaluation
    if (classifier_config.run_evaluation ?: false) {
        evaluation_result = EVALUATE_CLASSIFIER(
            processed_data_ch,
            db_type,
            classifier_config.n_jobs ?: 1
        )
        
        evaluation_ch = evaluation_result.evaluation
        predicted_taxonomy_ch = evaluation_result.predicted_taxonomy
    } else {
        evaluation_ch = Channel.empty()
        predicted_taxonomy_ch = Channel.empty()
    }

    emit:
    classifier = classifier_ch
    read_extraction_stats = processed_result.read_extraction_stats
    evaluation = evaluation_ch
    predicted_taxonomy = predicted_taxonomy_ch
    processed_sequences = processed_data_ch.map { _meta, seqs, tax -> tuple(meta, seqs) }
    processed_taxonomy = processed_data_ch.map { _meta, seqs, tax -> tuple(meta, tax) }
}

// Convenience wrapper workflows for backwards compatibility

workflow CREATE_SILVA_CLASSIFIER {
    take:
    meta
    silva_version
    silva_target
    filter_labels
    filter_min_lengths
    derep_mode
    f_primer
    r_primer
    run_evaluation

    main:
    db_config = [
        version: silva_version,
        target: silva_target,
        include_species_labels: false
    ]
    
    processing_config = [
        cull_sequences: true,
        filter_by_length: (filter_labels && filter_min_lengths),
        filter_labels: filter_labels,
        filter_min_lengths: filter_min_lengths,
        derep_mode: derep_mode,
        f_primer: f_primer,
        r_primer: r_primer
    ]
    
    classifier_config = [
        run_evaluation: run_evaluation ?: false,
        n_jobs: 1
    ]

    result = CREATE_REF_DB_CLASSIFIER(
        meta,
        'silva',
        db_config,
        processing_config,
        classifier_config
    )

    emit:
    classifier = result.classifier
    evaluation = result.evaluation
    ref_seqs   = result.processed_sequences.map { _meta, seqs -> seqs }
    read_extraction_stats = result.read_extraction_stats
}

workflow CREATE_UNITE_CLASSIFIER {
    take:
    meta
    unite_version
    unite_taxon_group
    unite_cluster_id
    include_singletons
    its_region
    derep_mode
    run_evaluation

    main:
    db_config = [
        version: unite_version,
        taxon_group: unite_taxon_group,
        cluster_id: unite_cluster_id,
        include_singletons: include_singletons ?: false
    ]
    
    processing_config = [
        remove_sh_identifiers: true,
        its_region: its_region,
        derep_mode: derep_mode
    ]
    
    classifier_config = [
        run_evaluation: run_evaluation ?: false,
        n_jobs: 1
    ]

    result = CREATE_REF_DB_CLASSIFIER(
        meta,
        'unite',
        db_config,
        processing_config,
        classifier_config
    )

    emit:
    classifier = result.classifier
    evaluation = result.evaluation
}

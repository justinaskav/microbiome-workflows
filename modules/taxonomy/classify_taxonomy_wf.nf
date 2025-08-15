// Workflow for Taxonomic Classification and Visualization using QIIME2

// Include visualization workflow
include { VISUALIZE_TAXONOMY } from './taxonomy_visualization_wf.nf'

// Default parameters specific to this workflow module
// These can be overridden in nextflow.config or main script
params.classifier_confidence = 0.7        // Confidence threshold for taxonomic classification
params.classifier_read_orientation = 'auto' // Read orientation for classification
params.classifier_cpus = 2                // Default CPUs for classification process

// Process to classify representative sequences using a pre-trained sklearn classifier
process CLASSIFY_SKLEARN_QIIME2 {
    tag "classify_taxonomy_${rep_seqs_qza.getBaseName()}"
    conda "${params.qiime2_conda_env}"
    cpus params.classifier_cpus ?: 4 // Use module-specific or default CPUs

    input:
    path rep_seqs_qza     // Representative sequences (QZA)
    path classifier_qza   // Pre-trained classifier (QZA)

    output:
    path("taxonomy.qza"), emit: classification_qza

    script:
    def n_jobs = task.cpus > 0 ? task.cpus : 1 // Use task.cpus for --p-n-jobs

    """
    qiime feature-classifier classify-sklearn \
        --i-reads "${rep_seqs_qza}" \
        --i-classifier "${classifier_qza}" \
        --p-confidence ${params.classifier_confidence} \
        --p-read-orientation "${params.classifier_read_orientation}" \
        --p-n-jobs ${n_jobs} \
        --o-classification taxonomy.qza \
        --verbose
    """
}

workflow CLASSIFY_TAXONOMY {
    take:
        rep_seqs_ch      // Channel of representative sequences (QZA)
        classifier_ch    // Channel of the pre-trained classifier (QZA)
        feature_table_ch // Channel of the feature table (QZA) for visualization
        metadata_ch      // Channel of metadata (optional, for visualization)

    main:
        // 1. Classify representative sequences
        CLASSIFY_SKLEARN_QIIME2 (
            rep_seqs_ch,
            classifier_ch
        )
        classification_results_ch = CLASSIFY_SKLEARN_QIIME2.out.classification_qza

        // 2. Visualize taxonomic assignments using the common visualization workflow
        taxonomy_viz_ch = VISUALIZE_TAXONOMY(
            classification_results_ch,
            feature_table_ch,
            "taxonomy_barplot",
            metadata_ch ? metadata_ch : null // Use first metadata file if available
        )
        taxonomy_summary_ch = taxonomy_viz_ch.barplot

    emit:
        classification = classification_results_ch   // The .qza file with taxonomic assignments
        classification_summary = taxonomy_summary_ch // The .qzv visualization
}


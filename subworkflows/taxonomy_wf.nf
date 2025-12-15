// Taxonomy Subworkflow
// Performs taxonomic classification with optional DADA2 taxonomy, comparison, and gram staining evaluation

nextflow.enable.dsl = 2

include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { EVALUATE_REFERENCE_TAXONOMY } from '../modules/reference_databases/evaluate_taxonomy_wf.nf'
include { COMPARE_TAXONOMY_ARTIFACTS } from '../modules/taxonomy/taxonomy_comparison_wf.nf'
include { WF_EVALUATE_GRAM_STAINING } from '../modules/taxonomy/evaluate_gram_staining_wf.nf'

workflow TAXONOMY {
    take:
    rep_seqs           // Filtered representative sequences
    feature_table      // Filtered feature table
    metadata_file      // Metadata file
    database_ch_list   // Map with classifier and refs (trainset, species)
    taxonomy_options   // Map with enable_dada2_taxonomy, enable_comparison

    main:

    taxonomy_config = [
        enable_dada2_taxonomy: taxonomy_options?.enable_dada2_taxonomy ?: params.taxonomy?.enable_dada2_taxonomy ?: false,
        enable_comparison: taxonomy_options?.enable_comparison ?: params.taxonomy?.enable_comparison ?: false,
    ]

    // Taxonomic Classification with feature-classifier
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        rep_seqs,
        database_ch_list.classifier,
        feature_table,
        metadata_file,
    )

    classified_taxa_ch = taxonomy_results_ch.classification
    classified_taxa_summary_ch = taxonomy_results_ch.classification_summary

    // Optional DADA2 taxonomy assignment
    if (taxonomy_config.enable_dada2_taxonomy) {
        dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
            rep_seqs,
            feature_table,
            database_ch_list.refs.trainset,
            database_ch_list.refs.species,
            metadata_file
        )

        if (taxonomy_config.enable_comparison) {
            // Evaluate taxonomic classifications using RESCRIPt
            taxonomy_eval_ch = EVALUATE_REFERENCE_TAXONOMY(
                classified_taxa_ch.mix(dada2_taxonomy_results_ch.qiime2_taxonomy),
                "taxonomy_evaluation",
                "",
            ).evaluation

            taxonomy_comparison_ch_results = COMPARE_TAXONOMY_ARTIFACTS(
                classified_taxa_ch,
                dada2_taxonomy_results_ch.qiime2_taxonomy,
                Channel.value("Feature-Classifier"),
                Channel.value("DADA2"),
            )

            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
            dada2_taxonomy_summary_ch = dada2_taxonomy_results_ch.taxonomy_barplot
            taxonomy_evaluation_ch = taxonomy_eval_ch
            taxonomy_comparison_ch = taxonomy_comparison_ch_results.comparison_results
        }
        else {
            dada2_taxonomy_ch = dada2_taxonomy_results_ch.qiime2_taxonomy
            dada2_taxonomy_summary_ch = dada2_taxonomy_results_ch.taxonomy_barplot
            taxonomy_evaluation_ch = Channel.empty()
            taxonomy_comparison_ch = Channel.empty()
        }
    }
    else {
        // Skip DADA2 taxonomy - create empty channels
        dada2_taxonomy_ch = Channel.empty()
        dada2_taxonomy_summary_ch = Channel.empty()
        taxonomy_evaluation_ch = Channel.empty()
        taxonomy_comparison_ch = Channel.empty()
    }

    // Gram Staining Evaluation
    WF_EVALUATE_GRAM_STAINING(
        classified_taxa_ch,
        feature_table,
        metadata_file,
    )

    emit:
    taxonomy = classified_taxa_ch
    taxonomy_summary = classified_taxa_summary_ch
    dada2_taxonomy = dada2_taxonomy_ch
    dada2_taxonomy_summary = dada2_taxonomy_summary_ch
    taxonomy_evaluation = taxonomy_evaluation_ch
    taxonomy_comparison = taxonomy_comparison_ch
    gram_staining_evaluation = WF_EVALUATE_GRAM_STAINING.out.gram_staining_eval_results
}

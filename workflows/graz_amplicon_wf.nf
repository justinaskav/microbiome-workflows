// Reusable Graz Lab Amplicon Processing Workflow
// All parameters should be set in project-specific nextflow.config files

nextflow.enable.dsl = 2

include { DOWNLOAD_SILVA } from '../modules/silva/download_silva_wf.nf'
include { PRECLASSIFY_SILVA } from '../modules/silva/preclassify_silva_wf.nf'
include { WF_QUALITY_CONTROL } from '../modules/quality_control/quality_control_wf.nf'
include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'
include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'
include { EVALUATE_TAXONOMY } from '../modules/silva/rescript_evaluate_taxonomy_wf.nf'
include { COMPARE_TAXONOMY_ARTIFACTS } from '../modules/taxonomy/taxonomy_comparison_wf.nf'

// Named workflow for Graz lab amplicon processing
workflow GRAZ_AMPLICON_PROCESSING {
    take:
    reads_pattern              // File pattern for paired-end reads
    silva_version             // SILVA database version
    silva_target              // SILVA target (e.g., 'SSURef_NR99')
    filter_labels             // Filter labels for SILVA
    filter_min_lengths        // Filter minimum lengths
    derep_mode                // Dereplication mode
    f_primer                  // Forward primer
    r_primer                  // Reverse primer
    dada2_trunc_len_f         // DADA2 forward truncation length
    dada2_trunc_len_r         // DADA2 reverse truncation length
    dada2_trim_left_f         // DADA2 forward trim left
    dada2_trim_left_r         // DADA2 reverse trim left
    dada2_n_reads_learn       // DADA2 reads for learning
    dada2_n_threads           // DADA2 threads

    main:
    // Create read pairs channel
    read_pairs_ch = Channel.fromFilePairs(reads_pattern)

    // Quality Control
    WF_QUALITY_CONTROL(read_pairs_ch)

    // Download and prepare SILVA classifier
    def silva_meta = [id: "silva-${silva_version}-${silva_target}"]
    
    silva_classifier_ch = PRECLASSIFY_SILVA(
        silva_meta,
        silva_version,
        silva_target,
        filter_labels,
        filter_min_lengths,
        derep_mode,
        f_primer,
        r_primer
    ).classifier

    // DADA2 denoising
    dada2_output_ch = WF_DENOISE_DADA2(
        read_pairs_ch,
        dada2_trunc_len_f,
        dada2_trunc_len_r,
        dada2_trim_left_f,
        dada2_trim_left_r,
        dada2_n_reads_learn,
        dada2_n_threads
    )

    // Taxonomic classification
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        dada2_output_ch.representative_sequences,
        silva_classifier_ch,
        dada2_output_ch.feature_table,
        Channel.empty() // No metadata filter for basic workflow
    )

    // Alternative DADA2 taxonomy
    dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
        dada2_output_ch.representative_sequences,
        dada2_output_ch.feature_table
    )

    // Phylogenetic analysis
    ALIGN_AND_TREE_QIIME2(dada2_output_ch.representative_sequences)

    // Taxonomy evaluation and comparison
    taxonomy_eval_ch = EVALUATE_TAXONOMY(
        taxonomy_results_ch.classification.map { it -> [it] }
            .mix(dada2_taxonomy_results_ch.qiime2_taxonomy.map { it -> [it] }),
        Channel.value("taxonomy_evaluation")
    )

    compare_script_ch = Channel.fromPath("../scripts/compare_taxonomies.py")
    taxonomy_comparison_ch = COMPARE_TAXONOMY_ARTIFACTS(
        taxonomy_results_ch.classification,
        dada2_taxonomy_results_ch.qiime2_taxonomy,
        Channel.value("Feature-Classifier"),
        Channel.value("DADA2"),
        compare_script_ch
    )

    emit:
    // Core outputs
    feature_table = dada2_output_ch.feature_table
    rep_seqs = dada2_output_ch.representative_sequences
    taxonomy = taxonomy_results_ch.classification
    taxonomy_summary = taxonomy_results_ch.classification_summary
    phylogeny_rooted = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
    phylogeny_unrooted = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
    
    // QC outputs
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    
    // Alternative taxonomy
    dada2_taxonomy = dada2_taxonomy_results_ch.qiime2_taxonomy
    
    // Evaluation outputs
    taxonomy_evaluation = taxonomy_eval_ch.taxonomy_stats
    taxonomy_comparison = taxonomy_comparison_ch.comparison_results
    
    // Classifier
    classifier = silva_classifier_ch
}
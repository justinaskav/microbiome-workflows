// Amplicon Post-Trim Workflow
// Takes trimmed reads (potentially merged from multiple batches) and runs DADA2 → downstream
// Used for the pre-DADA2 merge strategy

nextflow.enable.dsl = 2

// Include subworkflows
include { DADA2 } from '../subworkflows/dada2_wf.nf'
include { FILTERING } from '../subworkflows/filtering_wf.nf'
include { TAXONOMY } from '../subworkflows/taxonomy_wf.nf'
include { PICRUST2 } from '../subworkflows/picrust2_wf.nf'
include { PHYLOGENY } from '../subworkflows/phylogeny_wf.nf'
include { DIVERSITY } from '../subworkflows/diversity_wf.nf'
include { WF_FILTER_SAMPLES_PRE_DADA2 } from '../modules/data_processing/filter_samples_pre_dada2_wf.nf'

workflow AMPLICON_POST_TRIM {
    take:
    trimmed_reads_ch      // Channel of [sample_id, [r1, r2]] - already trimmed
    metadata_file         // Merged metadata file
    dada2_options
    decontam_options
    taxonomy_options
    database_ch_list
    other_options

    main:

    // Step 0.5: Pre-DADA2 sample filtering (optional)
    sequences_for_dada2 = trimmed_reads_ch
    if (other_options?.filter_samples_by_metadata) {
        pre_filter_results = WF_FILTER_SAMPLES_PRE_DADA2(
            sequences_for_dada2,
            metadata_file,
        )
        sequences_for_dada2 = pre_filter_results.filtered_reads
    }

    // Step 1: DADA2 Denoising on trimmed/merged reads
    dada2_results = DADA2(sequences_for_dada2, dada2_options)

    // Step 2: Filtering (metadata filtering + decontamination)
    filtering_results = FILTERING(
        dada2_results.feature_table,
        dada2_results.rep_seqs,
        metadata_file,
        decontam_options,
        database_ch_list.classifier,
    )

    // Step 3: Taxonomy Classification
    taxonomy_results = TAXONOMY(
        filtering_results.filtered_rep_seqs,
        filtering_results.filtered_feature_table,
        metadata_file,
        database_ch_list,
        taxonomy_options
    )

    // Step 4: PICRUSt2 Functional Prediction (skipped for ITS)
    picrust2_ko_metagenome_ch = Channel.empty()
    picrust2_ec_metagenome_ch = Channel.empty()
    picrust2_pathway_abundance_ch = Channel.empty()
    picrust2_ko_summary_ch = Channel.empty()
    picrust2_ec_summary_ch = Channel.empty()
    picrust2_pathway_summary_ch = Channel.empty()
    picrust2_pathway_jaccard_emperor_ch = Channel.empty()
    picrust2_pathway_bray_curtis_emperor_ch = Channel.empty()

    if (!other_options?.skip_picrust2) {
        picrust2_options = [
            threads: other_options?.picrust2_threads,
            diversity_sampling_depth: other_options?.picrust2_diversity_sampling_depth,
            skip_pathway_diversity: other_options?.skip_picrust2_diversity ?: false,
            auto_depth_percentile: other_options?.auto_depth_percentile ?: 0.0,
        ]

        picrust2_results = PICRUST2(
            filtering_results.filtered_feature_table,
            filtering_results.filtered_rep_seqs,
            metadata_file,
            picrust2_options,
        )

        picrust2_ko_metagenome_ch = picrust2_results.ko_metagenome
        picrust2_ec_metagenome_ch = picrust2_results.ec_metagenome
        picrust2_pathway_abundance_ch = picrust2_results.pathway_abundance
        picrust2_ko_summary_ch = picrust2_results.ko_summary
        picrust2_ec_summary_ch = picrust2_results.ec_summary
        picrust2_pathway_summary_ch = picrust2_results.pathway_summary
        picrust2_pathway_jaccard_emperor_ch = picrust2_results.pathway_jaccard_emperor
        picrust2_pathway_bray_curtis_emperor_ch = picrust2_results.pathway_bray_curtis_emperor
    }

    // Step 5: Phylogeny
    phylogeny_results = PHYLOGENY(filtering_results.filtered_rep_seqs)

    // Step 6: Diversity Analysis (optional)
    diversity_options = [
        sampling_depth: other_options?.diversity_sampling_depth,
        amplicon_type: other_options?.amplicon_type,
        database: other_options?.database,
        skip: other_options?.skip_diversity ?: false,
        auto_depth_percentile: other_options?.auto_depth_percentile ?: 0.0,
    ]

    diversity_results = DIVERSITY(
        filtering_results.filtered_feature_table,
        phylogeny_results.rooted_tree,
        metadata_file,
        diversity_options
    )

    emit:
    // DADA2 results
    dada2_denoising_stats = dada2_results.denoising_stats
    dada2_ft_summary = dada2_results.feature_table_summary
    dada2_rep_seqs_summary = dada2_results.rep_seqs_summary
    dada2_denoising_stats_summary = dada2_results.denoising_stats_summary
    dada2_base_transition_stats = dada2_results.base_transition_stats

    // Filtering/Decontam results
    decontam_scores = filtering_results.decontam_scores
    decontam_scores_viz = filtering_results.decontam_scores_viz
    contaminant_table = filtering_results.contaminant_table
    contaminant_seqs = filtering_results.contaminant_seqs
    contaminant_taxonomy = filtering_results.contaminant_taxonomy

    // Taxonomy results
    taxonomy = taxonomy_results.taxonomy
    taxonomy_summary = taxonomy_results.taxonomy_summary
    dada2_taxonomy = taxonomy_results.dada2_taxonomy
    dada2_taxonomy_summary = taxonomy_results.dada2_taxonomy_summary
    taxonomy_evaluation = taxonomy_results.taxonomy_evaluation
    taxonomy_comparison = taxonomy_results.taxonomy_comparison
    gram_staining_evaluation = taxonomy_results.gram_staining_evaluation

    // PICRUSt2 results
    picrust2_ko_metagenome = picrust2_ko_metagenome_ch
    picrust2_ec_metagenome = picrust2_ec_metagenome_ch
    picrust2_pathway_abundance = picrust2_pathway_abundance_ch
    picrust2_ko_summary = picrust2_ko_summary_ch
    picrust2_ec_summary = picrust2_ec_summary_ch
    picrust2_pathway_summary = picrust2_pathway_summary_ch
    picrust2_pathway_jaccard_emperor = picrust2_pathway_jaccard_emperor_ch
    picrust2_pathway_bray_curtis_emperor = picrust2_pathway_bray_curtis_emperor_ch

    // Phylogeny results
    phylogeny_rooted = phylogeny_results.rooted_tree
    phylogeny_unrooted = phylogeny_results.unrooted_tree

    // Diversity results
    rarefied_table = diversity_results.rarefied_table
    jaccard_emperor = diversity_results.jaccard_emperor
    bray_curtis_emperor = diversity_results.bray_curtis_emperor

    // Final outputs
    feature_table = filtering_results.filtered_feature_table
    rep_seqs = filtering_results.filtered_rep_seqs
}

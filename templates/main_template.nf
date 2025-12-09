#!/usr/bin/env nextflow
// Template main.nf demonstrating both merge strategies
// Strategy 1: Post-DADA2 merge (merge feature tables after individual DADA2 runs)
// Strategy 2: Pre-DADA2 merge (merge trimmed reads, run single DADA2)

nextflow.enable.dsl = 2

// === IMPORTS ===

// Reference database workflows
include { DOWNLOAD_DADA2_SILVA_REFERENCES } from './microbiome-workflows/modules/reference_databases/download_reference_db_wf.nf'
include { CREATE_SILVA_CLASSIFIER } from './microbiome-workflows/modules/reference_databases/create_classifier_wf.nf'

// Full amplicon workflow (for post-DADA2 merge strategy)
include { AMPLICON_WORKFLOW as BATCH_1 ; AMPLICON_WORKFLOW as BATCH_2 } from './microbiome-workflows/workflows/amplicon_wf.nf'

// Trimming-only workflow (for pre-DADA2 merge strategy)
include { AMPLICON_TRIMMING_ONLY as TRIM_BATCH_1 ; AMPLICON_TRIMMING_ONLY as TRIM_BATCH_2 } from './microbiome-workflows/workflows/amplicon_trimming_only_wf.nf'

// Post-trim workflow (for pre-DADA2 merge strategy)
include { AMPLICON_POST_TRIM } from './microbiome-workflows/workflows/amplicon_post_trim_wf.nf'

// Merge workflows
include { WF_MERGE_TRIMMED_READS ; WF_MERGE_BATCH_METADATA } from './microbiome-workflows/modules/data_processing/merge_trimmed_reads_wf.nf'
include { MERGE_FEATURE_TABLES ; MERGE_SEQUENCES ; MERGE_TAXONOMIES } from './microbiome-workflows/modules/data_processing/merge_qiime2_wf.nf'

// Post-merge analysis
include { WF_CREATE_SUMMARIES } from './microbiome-workflows/modules/visualization/summary_visualization_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from './microbiome-workflows/modules/phylogeny/phylogeny_wf.nf'
include { WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC } from './microbiome-workflows/modules/diversity/diversity_core_metrics_wf.nf'


// === MAIN WORKFLOW ===

workflow {
    main:

    // Download DADA2 references (shared by all strategies)
    silva_refs = DOWNLOAD_DADA2_SILVA_REFERENCES()

    // Create classifier (adjust primers per your region)
    classifier_ch = CREATE_SILVA_CLASSIFIER(
        [id: "silva-${params.silva.version}-V3V4"],
        params.silva.version ?: '138.2',
        params.silva.target ?: 'SSURef_NR99',
        params.silva.filter_labels ?: ['Archaea', 'Bacteria', 'Eukaryota'],
        params.silva.filter_min_lengths ?: [900, 1200, 1400],
        params.silva.derep_mode ?: 'uniq',
        params.batches[0].primers.forward,
        params.batches[0].primers.reverse,
        false,
    ).classifier

    database_refs = [
        classifier: classifier_ch,
        refs: [trainset: silva_refs.trainset, species: silva_refs.species]
    ]

    // =========================================================================
    // STRATEGY SELECTION
    // Set params.merge_strategy = "pre_dada2" or "post_dada2" in config
    // =========================================================================

    if (params.merge_strategy == "pre_dada2") {
        // =====================================================================
        // STRATEGY 2: PRE-DADA2 MERGE
        // Trim each batch separately, merge reads, run single DADA2
        // Better for: consistent ASV calling across batches, shared error model
        // =====================================================================

        // Step 1: Run trimming only on each batch
        trim1_results = TRIM_BATCH_1(
            params.batches[0].name,
            "${projectDir}/${params.batches[0].metadata_path}",
            params.batches[0].input,
            params.batches[0].primers,
            params.batches[0].trimming,
        )

        trim2_results = TRIM_BATCH_2(
            params.batches[1].name,
            "${projectDir}/${params.batches[1].metadata_path}",
            params.batches[1].input,
            params.batches[1].primers,
            params.batches[1].trimming,
        )

        // Step 2: Merge trimmed reads from all batches
        merged_reads_ch = WF_MERGE_TRIMMED_READS(
            trim1_results.trimmed_reads.mix(trim2_results.trimmed_reads)
        ).merged_reads

        // Step 3: Merge metadata files
        merged_metadata = WF_MERGE_BATCH_METADATA(
            trim1_results.metadata_file.mix(trim2_results.metadata_file)
        ).merged_metadata

        // Step 4: Run DADA2 and downstream on merged data
        // Use DADA2 params appropriate for merged data (typically same as batch 1)
        merged_results = AMPLICON_POST_TRIM(
            merged_reads_ch,
            merged_metadata,
            params.batches[0].dada2,  // or create unified dada2 params
            params.decontam,
            params.taxonomy,
            database_refs,
            [diversity_sampling_depth: params.diversity_sampling_depth],
        )

        // Assign outputs for the output block
        feature_table = merged_results.feature_table
        rep_seqs = merged_results.rep_seqs
        taxonomy = merged_results.taxonomy
        rooted_tree = merged_results.phylogeny_rooted

    } else {
        // =====================================================================
        // STRATEGY 1: POST-DADA2 MERGE (Default)
        // Run full pipeline on each batch, merge feature tables at the end
        // Better for: batches with different quality profiles, existing workflows
        // =====================================================================

        // Process each batch through full pipeline
        batch1_results = BATCH_1(
            params.batches[0].name,
            "${projectDir}/${params.batches[0].metadata_path}",
            params.batches[0].input,
            params.batches[0].primers,
            params.batches[0].trimming,
            params.batches[0].dada2,
            params.batches[0].decontam,
            params.taxonomy,
            database_refs,
            params.batches[0].other ?: [:],
        )

        batch2_results = BATCH_2(
            params.batches[1].name,
            "${projectDir}/${params.batches[1].metadata_path}",
            params.batches[1].input,
            params.batches[1].primers,
            params.batches[1].trimming,
            params.batches[1].dada2,
            params.batches[1].decontam,
            params.taxonomy,
            database_refs,
            params.batches[1].other ?: [:],
        )

        // Collect outputs from all batches
        feature_tables_ch = batch1_results.feature_table
            .mix(batch2_results.feature_table)
            .collect()
        rep_seqs_ch = batch1_results.rep_seqs
            .mix(batch2_results.rep_seqs)
            .collect()
        taxonomy_ch = batch1_results.taxonomy
            .mix(batch2_results.taxonomy)
            .collect()

        // Merge results
        merged_tables = MERGE_FEATURE_TABLES(feature_tables_ch)
        merged_seqs = MERGE_SEQUENCES(rep_seqs_ch)
        merged_taxa = MERGE_TAXONOMIES(taxonomy_ch)

        // Assign outputs
        feature_table = merged_tables.merged_table
        rep_seqs = merged_seqs.merged_seqs
        taxonomy = merged_taxa.merged_taxonomy

        // Generate phylogeny from merged sequences
        phylogeny = ALIGN_AND_TREE_QIIME2(merged_seqs.merged_seqs)
        rooted_tree = phylogeny.rooted_tree_qza
    }

    // Common post-merge analysis
    metadata_file = Channel.fromPath(params.metadata_file)

    summaries = WF_CREATE_SUMMARIES(
        feature_table,
        rep_seqs,
        metadata_file,
    )

    diversity = WF_DIVERSITY_CORE_METRICS_PHYLOGENETIC(
        feature_table,
        rooted_tree,
        metadata_file,
        params.diversity_sampling_depth ?: 10000,
        [id: "merged"],
    )

    // === WORKFLOW OUTPUTS ===
    // These feed into the output {} block below

    merged_feature_table = feature_table
    merged_rep_seqs = rep_seqs
    merged_taxonomy = taxonomy
    merged_rooted_tree = rooted_tree
    merged_summaries = summaries
    merged_diversity = diversity
}

// === OUTPUT BLOCK ===
// Defines where outputs are published

output {
    // Merged outputs
    merged_feature_table {
        path 'merged/feature_table'
    }
    merged_rep_seqs {
        path 'merged/rep_seqs'
    }
    merged_taxonomy {
        path 'merged/taxonomy'
    }
    merged_rooted_tree {
        path 'merged/phylogeny'
    }

    // Summaries
    'merged_summaries.feature_table_summary' {
        path 'merged/summaries'
    }
    'merged_summaries.rep_seqs_summary' {
        path 'merged/summaries'
    }

    // Diversity
    'merged_diversity.jaccard_emperor' {
        path 'merged/diversity'
    }
    'merged_diversity.bray_curtis_emperor' {
        path 'merged/diversity'
    }
    'merged_diversity.rarefied_table' {
        path 'merged/diversity'
    }
}

// Amplicon Workflow - Main Orchestrator
// Chains all subworkflows for complete amplicon analysis pipeline

nextflow.enable.dsl = 2

// Include subworkflows
include { INPUT_PROCESSING } from '../subworkflows/input_processing_wf.nf'
include { QC } from '../subworkflows/qc_wf.nf'
include { TRIMMING } from '../subworkflows/trimming_wf.nf'
include { DADA2 } from '../subworkflows/dada2_wf.nf'
include { DEBLUR } from '../subworkflows/deblur_wf.nf'
include { FILTERING } from '../subworkflows/filtering_wf.nf'
include { TAXONOMY } from '../subworkflows/taxonomy_wf.nf'
include { PICRUST2 } from '../subworkflows/picrust2_wf.nf'
include { PHYLOGENY } from '../subworkflows/phylogeny_wf.nf'
include { DIVERSITY } from '../subworkflows/diversity_wf.nf'
include { WF_FILTER_SAMPLES_PRE_DADA2 } from '../modules/data_processing/filter_samples_pre_dada2_wf.nf'

// Main amplicon workflow
workflow AMPLICON_WORKFLOW {
    take:
    name
    metadata_path
    input_options
    primers
    trimming_options
    dada2_options
    decontam_options
    taxonomy_options
    database_ch_list
    other_options

    main:

    // Step 1: Input Processing
    input_processing_results = INPUT_PROCESSING(input_options)

    // Use dynamic metadata if available, otherwise use static path
    metadata_to_use = input_processing_results.metadata_file.ifEmpty(
        metadata_path ? file(metadata_path) : Channel.empty()
    )

    // Step 2: Quality Control on original reads
    qc_results = QC(input_processing_results.read_pairs_ch)

    // Step 3: Trimming (optional)
    trimming_config = [
        enabled: trimming_options.enabled ?: false,
        use_sampleseq: trimming_options.use_sampleseq ?: false,
        skip_polyg_trim: trimming_options.skip_polyg_trim ?: false,
    ]

    if (trimming_config.enabled) {
        trimming_results = TRIMMING(
            input_processing_results.read_pairs_ch,
            input_processing_results.barcode_primer_info_ch,
            primers,
            trimming_options
        )
        sequences_for_dada2 = trimming_results.trimmed_reads
        fastqc_after_cutadapt = trimming_results.fastqc_results_after_trim
        multiqc_after_cutadapt = trimming_results.multiqc_report_after_trim
    } else {
        sequences_for_dada2 = input_processing_results.read_pairs_ch
        fastqc_after_cutadapt = Channel.empty()
        multiqc_after_cutadapt = Channel.empty()
    }

    // Step 3.5: Pre-DADA2 sample filtering (optional)
    if (other_options?.filter_samples_by_metadata) {
        pre_filter_results = WF_FILTER_SAMPLES_PRE_DADA2(
            sequences_for_dada2,
            metadata_to_use,
        )
        sequences_for_dada2 = pre_filter_results.filtered_reads
    }

    // Step 4: Denoising (DADA2 or Deblur based on config)
    def denoising_method = other_options?.denoising_method ?: 'dada2'

    if (denoising_method == 'deblur') {
        denoising_results = DEBLUR(
            sequences_for_dada2,
            other_options.deblur,
            database_ch_list.deblur_ref_seqs ?: Channel.empty(),
        )
        denoising_base_transition_ch = Channel.empty()
    } else {
        denoising_results = DADA2(sequences_for_dada2, dada2_options)
        denoising_base_transition_ch = denoising_results.base_transition_stats
    }

    // Step 5: Filtering (metadata filtering + decontamination)
    // Pass classifier for contaminant taxonomy classification
    filtering_results = FILTERING(
        denoising_results.feature_table,
        denoising_results.rep_seqs,
        metadata_to_use,
        decontam_options,
        database_ch_list.classifier
    )

    // Step 6: Taxonomy Classification
    taxonomy_results = TAXONOMY(
        filtering_results.filtered_rep_seqs,
        filtering_results.filtered_feature_table,
        metadata_to_use,
        database_ch_list,
        taxonomy_options
    )

    // Step 7: PICRUSt2 Functional Prediction (skipped for ITS)
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
            metadata_to_use,
            picrust2_options
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

    // Step 8: Phylogeny
    phylogeny_results = PHYLOGENY(filtering_results.filtered_rep_seqs)

    // Step 9: Diversity Analysis (optional)
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
        metadata_to_use,
        diversity_options
    )

    emit:
    // QC results
    fastqc_results = qc_results.fastqc_results
    multiqc_report = qc_results.multiqc_report
    fastqc_results_after_cutadapt = fastqc_after_cutadapt
    multiqc_report_after_cutadapt = multiqc_after_cutadapt

    // Denoising results (DADA2 or Deblur)
    dada2_denoising_stats = denoising_results.denoising_stats
    dada2_ft_summary = denoising_results.feature_table_summary
    dada2_rep_seqs_summary = denoising_results.rep_seqs_summary
    dada2_denoising_stats_summary = denoising_results.denoising_stats_summary
    dada2_base_transition_stats_summary = denoising_base_transition_ch

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
    phylogeny_unrooted = phylogeny_results.unrooted_tree
    phylogeny_rooted = phylogeny_results.rooted_tree

    // Diversity results
    rarefied_table = diversity_results.rarefied_table
    jaccard_emperor = diversity_results.jaccard_emperor
    bray_curtis_emperor = diversity_results.bray_curtis_emperor

    // Final outputs
    feature_table = filtering_results.filtered_feature_table
    rep_seqs = filtering_results.filtered_rep_seqs
}

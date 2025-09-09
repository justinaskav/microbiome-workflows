// Reusable Amplicon Processing Workflow
// All parameters should be set in project-specific nextflow.config files

nextflow.enable.dsl = 2

// Include required workflows and modules
include { PRECLASSIFY_SILVA } from '../modules/silva/preclassify_silva_wf.nf'
include { WF_QUALITY_CONTROL; WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_CUTADAPT } from '../modules/quality_control/quality_control_wf.nf'
include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_PARSE_NOVOGENE_SAMPLESEQ } from '../modules/novogene/sampleseq_wf.nf'
include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'
include { WF_FILTER_BY_METADATA } from '../modules/data_processing/filter_by_metadata_wf.nf'
include { WF_DECONTAM_FREQUENCY; WF_DECONTAM_PREVALENCE } from '../modules/decontam/decontam_wf.nf'
include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'
include { EVALUATE_TAXONOMY } from '../modules/silva/rescript_evaluate_taxonomy_wf.nf'
include { COMPARE_TAXONOMY_ARTIFACTS } from '../modules/taxonomy/taxonomy_comparison_wf.nf'
include { WF_UNZIP_DATA } from '../modules/data_processing/unzip_data_wf.nf'

// Named workflow that can be called from other workflows
workflow AMPLICON_WORKFLOW {
    take:
    data_dir                    // Data directory path
    zip_filename               // Zip filename (without .zip extension), null for unzipped data
    metadata_file             // Metadata file path
    batch_dada2_trunc_len_f   // DADA2 forward truncation length for this batch
    batch_dada2_trunc_len_r   // DADA2 reverse truncation length for this batch
    batch_dada2_trim_left_f // DADA2 forward trim left for this batch
    batch_dada2_trim_left_r // DADA2 reverse trim left for this batch
    batch_samples_with_additional_barcodes // Samples with additional barcodes for this batch
    silva_classifier_ch       // Pre-built SILVA classifier channel
    use_frequency_decontam    // Use frequency-based decontamination (optional, defaults to true)
    use_prevalence_decontam   // Use prevalence-based decontamination (optional, defaults to false)
    primer_f                  // Forward primer (for unzipped data), null to use SampleSeq from zip
    primer_r                  // Reverse primer (for unzipped data), null to use SampleSeq from zip
    skip_polyg_trim
    skip_dada2_taxonomy       // Skip DADA2 alternative taxonomy (optional, defaults to false)
    sample_prefix_to_remove   // Optional prefix to remove from sample IDs (e.g., 'FDMP20H000840-1a_L1_')

    main:
    // Set up metadata channel
    metadata_ch = Channel.fromPath(metadata_file)

    // Set decontamination parameters (use passed values or defaults)
    def use_freq_decontam = use_frequency_decontam != null ? use_frequency_decontam : params.use_decontam_frequency
    def use_prev_decontam = use_prevalence_decontam != null ? use_prevalence_decontam : params.use_decontam_prevalence

    // Handle both zipped and unzipped data
    if (zip_filename != null) {
        // Data is zipped - use the original unzip workflow
        def path_to_zip = file("${data_dir}/${zip_filename}.zip")
        unzip_results = WF_UNZIP_DATA(path_to_zip, zip_filename)
        fastq_files_ch = unzip_results.fastq_files
        sampleseq_file_ch = unzip_results.sampleseq_file

        // Process fastq files into read pairs for quality control
        fastq_files_ch
        .flatMap { list_of_files -> // Use flatMap to process each file in the list
            list_of_files.collect { file -> // 'file' here is a Path object
                def file_name = file.getName() // Use .getName() for Path objects
                def sample_name = file_name.replaceAll(/_(?:S[0-9]+_)?L001_R[12]_001\.fastq\.gz$/, "") // Handle optional S number
                return tuple(sample_name, file) // Emit a tuple: [sample_name, file_path]
            }
        }
        .groupTuple() // Group by sample name
        .set { read_pairs_ch } // Set the output channel
    } else {
        // Data is already in standard format - directly read fastq files
        // Look for standard paired-end fastq files
        read_pairs_ch = Channel.fromFilePairs("${data_dir}/*_L001_{R1,R2}_001.fastq.gz", size: 2) 

        // Trim _S part from sample names if present
        read_pairs_ch = read_pairs_ch.map { sample_name, files ->
            def trimmed_sample_name = sample_name.replaceAll(/_S\d+$/, "") // Remove any _S<number> suffix
            return tuple(trimmed_sample_name, files)
        }

        // Create empty sampleseq file channel since no SampleSeq_info.xls is needed
        sampleseq_file_ch = Channel.empty()
    }

    // Run regular quality control on original reads for comparison
    WF_QUALITY_CONTROL(read_pairs_ch)

    // Parse SampleSeq info if available (for Novogene data)
    if (zip_filename != null) {
        barcode_primer_info_ch = WF_PARSE_NOVOGENE_SAMPLESEQ(
            sampleseq_file_ch,
            batch_samples_with_additional_barcodes
        ).barcode_primer_info
    } else {
        barcode_primer_info_ch = Channel.empty()
    }

    // Run cutadapt workflow
    // Both zipped data (batch1/2) and standard format data (batch3) need cutadapt trimming
    trimming_results = WF_CUTADAPT(
        read_pairs_ch,
        barcode_primer_info_ch,
        true,  // Always enable cutadapt
        primer_f,
        primer_r,
        skip_polyg_trim
    )
    
    // Process exported fastqs for QC after cutadapt
    trimming_results.exported_fastqs
        .map { sample_name, fastq1, fastq2 -> 
            return tuple(sample_name, [fastq1, fastq2])
        }
        .set { exported_fastqs_pairs_ch }

    // Use exported fastqs for downstream analysis (DADA2 needs raw FASTQ files)
    trimmed_sequences_for_dada2 = exported_fastqs_pairs_ch

    // Run QC on the exported fastqs after cutadapt
    WF_QUALITY_CONTROL_AFTER_CUTADAPT(exported_fastqs_pairs_ch)

    // Use the passed SILVA classifier (created externally)

    dada2_output_ch = WF_DENOISE_DADA2(
        trimmed_sequences_for_dada2,
        batch_dada2_trunc_len_f,
        batch_dada2_trunc_len_r,
        batch_dada2_trim_left_f,
        batch_dada2_trim_left_r,
        params.dada2_n_reads_learn,
        params.dada2_n_threads,
        sample_prefix_to_remove
    )

    dada2_feature_table = dada2_output_ch.feature_table
    dada2_rep_seqs = dada2_output_ch.representative_sequences
    dada2_denoising_stats = dada2_output_ch.denoising_stats

    // Filter DADA2 outputs to only include samples present in metadata
    filtering_results_ch = WF_FILTER_BY_METADATA(
        dada2_feature_table,
        dada2_rep_seqs,
        metadata_ch
    )

    // Run decontamination based on enabled method
    if (use_freq_decontam) {
        decontam_results_ch = WF_DECONTAM_FREQUENCY(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_ch,
            params.decontam_freq_concentration_column,
            params.decontam_threshold
        )

        // Use decontaminated results for downstream analysis
        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores_freq = decontam_results_ch.decontam_scores
        decontam_scores_freq_viz = decontam_results_ch.decontam_scores_viz
        decontam_scores_prev = Channel.empty()
        decontam_scores_prev_viz = Channel.empty()
    } else if (use_prev_decontam) {
        decontam_results_ch = WF_DECONTAM_PREVALENCE(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_ch,
            params.decontam_prev_control_column,
            params.decontam_prev_control_indicator,
            params.decontam_threshold
        )

        // Use decontaminated results for downstream analysis
        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
        decontam_scores_prev = decontam_results_ch.decontam_scores
        decontam_scores_prev_viz = decontam_results_ch.decontam_scores_viz
        decontam_scores_freq = Channel.empty()
        decontam_scores_freq_viz = Channel.empty()
    } else {
        // Use metadata-filtered results directly for downstream analysis
        filtered_feature_table = filtering_results_ch.filtered_feature_table
        filtered_rep_seqs = filtering_results_ch.filtered_rep_seqs
        decontam_scores_freq = Channel.empty()
        decontam_scores_freq_viz = Channel.empty()
        decontam_scores_prev = Channel.empty()
        decontam_scores_prev_viz = Channel.empty()
    }

    // 4a. Taxonomic Classification with feature-classifier 
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        filtered_rep_seqs,
        silva_classifier_ch,
        filtered_feature_table,
        metadata_ch
    )

    classified_taxa_ch = taxonomy_results_ch.classification
    classified_taxa_summary_ch = taxonomy_results_ch.classification_summary

    // 4b. Alternative taxonomy classification using DADA2 (conditional)
    def run_dada2_taxonomy = skip_dada2_taxonomy != null ? !skip_dada2_taxonomy : true
    
    if (run_dada2_taxonomy) {
        dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
            filtered_rep_seqs,
            filtered_feature_table,
        )

        // 4c. Evaluate taxonomic classifications using RESCRIPt
        taxonomy_eval_ch = EVALUATE_TAXONOMY(
            classified_taxa_ch.map { it -> [it] }  // Convert to collection for multi-taxonomy support
                .mix(dada2_taxonomy_results_ch.qiime2_taxonomy.map { it -> [it] }),  // Add dada2 taxonomy
            Channel.value("taxonomy_evaluation") // Output filename
        )

        // 4d. Custom Pairwise Taxonomy Comparison
        compare_script_ch = Channel.fromPath("../scripts/compare_taxonomies.py")

        taxonomy_comparison_ch = COMPARE_TAXONOMY_ARTIFACTS(
            classified_taxa_ch,
            dada2_taxonomy_results_ch.qiime2_taxonomy,
            Channel.value("Feature-Classifier"),
            Channel.value("DADA2"),
            compare_script_ch
        )
    } else {
        // Skip DADA2 taxonomy - create empty channels
        dada2_taxonomy_results_ch = [qiime2_taxonomy: Channel.empty()]
        taxonomy_eval_ch = [taxonomy_stats: Channel.empty()]
        taxonomy_comparison_ch = [comparison_results: Channel.empty()]
    }

    // 5. Phylogenetic Analysis (Optional)
    ALIGN_AND_TREE_QIIME2(filtered_rep_seqs)

    emit:
    // Core outputs needed for merging
    feature_table = filtered_feature_table
    rep_seqs = filtered_rep_seqs
    taxonomy = classified_taxa_ch
    taxonomy_summary = classified_taxa_summary_ch

    // Additional outputs
    phylogeny_rooted = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
    phylogeny_unrooted = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
    dada2_taxonomy = dada2_taxonomy_results_ch.qiime2_taxonomy
    
    // QC outputs
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    fastqc_results_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.fastqc_results
    multiqc_report_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.multiqc_report
    
    // Denoising outputs
    dada2_denoising_stats = dada2_denoising_stats
    dada2_ft_summary = WF_DENOISE_DADA2.out.feature_table_summary
    dada2_rep_seqs_summary = WF_DENOISE_DADA2.out.rep_seqs_summary
    dada2_denoising_stats_summary = WF_DENOISE_DADA2.out.denoising_stats_summary

    // Decontamination outputs
    decontam_scores_freq = decontam_scores_freq
    decontam_scores_freq_viz = decontam_scores_freq_viz
    decontam_scores_prev = decontam_scores_prev
    decontam_scores_prev_viz = decontam_scores_prev_viz

    // Taxonomy comparison
    taxonomy_evaluation = taxonomy_eval_ch.taxonomy_stats
    taxonomy_comparison = taxonomy_comparison_ch.comparison_results
    
    // Classifier
    classifier = silva_classifier_ch
}


// Standardized QIIME2 Amplicon Processing Workflow
// Extracted from microbiome monorepo for reusability

nextflow.enable.dsl = 2

// Include required modules
include { WF_QUALITY_CONTROL; WF_QUALITY_CONTROL as WF_QUALITY_CONTROL_AFTER_CUTADAPT } from '../modules/quality_control/quality_control_wf.nf'
include { WF_CUTADAPT } from '../modules/sequencing/cutadapt_wf.nf'
include { WF_DENOISE_DADA2 } from '../modules/dada2/denoising_dada2_wf.nf'
include { WF_FILTER_BY_METADATA } from '../modules/data_processing/filter_by_metadata_wf.nf'
include { WF_DECONTAM_FREQUENCY; WF_DECONTAM_PREVALENCE } from '../modules/decontam/decontam_wf.nf'
include { CLASSIFY_TAXONOMY } from '../modules/taxonomy/classify_taxonomy_wf.nf'
include { WF_DADA2_TAXONOMY } from '../modules/dada2/dada2_assign_taxonomy_wf.nf'
include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'
include { WF_UNZIP_DATA } from '../modules/data_processing/unzip_data_wf.nf'
include { PRECLASSIFY_SILVA } from '../modules/silva/preclassify_silva_wf.nf'

// Main reusable workflow for QIIME2 amplicon processing
workflow QIIME2_AMPLICON_PROCESSING {
    take:
    data_dir                    // Data directory path
    zip_filename               // Zip filename (without .zip extension), null for unzipped data
    metadata_file             // Metadata file path
    dada2_trunc_len_f         // DADA2 forward truncation length
    dada2_trunc_len_r         // DADA2 reverse truncation length
    dada2_trim_left_f         // DADA2 forward trim left
    dada2_trim_left_r         // DADA2 reverse trim left
    samples_with_additional_barcodes // Samples with additional barcodes
    use_frequency_decontam    // Use frequency-based decontamination
    use_prevalence_decontam   // Use prevalence-based decontamination
    primer_f                  // Forward primer (optional)
    primer_r                  // Reverse primer (optional)
    silva_version            // SILVA database version (default: '138.2')
    silva_target             // SILVA target (default: 'SSURef_NR99')
    f_primer                 // Forward primer for classifier
    r_primer                 // Reverse primer for classifier

    main:
    // Set up metadata channel
    metadata_ch = Channel.fromPath(metadata_file)

    // Create SILVA classifier if not provided
    def silva_meta = [id: "silva-${silva_version}-${silva_target}"]
    silva_classifier_ch = PRECLASSIFY_SILVA(
        silva_meta,
        silva_version ?: '138.2',
        silva_target ?: 'SSURef_NR99',
        ['Archaea', 'Bacteria', 'Eukaryota'],
        [900, 1200, 1400],
        'uniq',
        f_primer,
        r_primer
    ).classifier

    // Handle both zipped and unzipped data
    if (zip_filename != null) {
        // Data is zipped - use the unzip workflow
        def path_to_zip = file("${data_dir}/${zip_filename}.zip")
        unzip_results = WF_UNZIP_DATA(path_to_zip, zip_filename)
        fastq_files_ch = unzip_results.fastq_files
        sampleseq_file_ch = unzip_results.sampleseq_file

        // Process fastq files into read pairs
        fastq_files_ch
        .flatMap { list_of_files ->
            list_of_files.collect { file ->
                def file_name = file.getName()
                def sample_name = file_name.replaceAll(/_(?:S[0-9]+_)?L001_R[12]_001\.fastq\.gz$/, "")
                return tuple(sample_name, file)
            }
        }
        .groupTuple()
        .set { read_pairs_ch }
    } else {
        // Data is already in standard format
        read_pairs_ch = Channel.fromFilePairs("${data_dir}/*_L001_{R1,R2}_001.fastq.gz", size: 2)
        read_pairs_ch = read_pairs_ch.map { sample_name, files ->
            def trimmed_sample_name = sample_name.replaceAll(/_S\d+$/, "")
            return tuple(trimmed_sample_name, files)
        }
        sampleseq_file_ch = Channel.empty()
    }

    // Run quality control on original reads
    WF_QUALITY_CONTROL(read_pairs_ch)

    // Run cutadapt workflow
    trimming_results = WF_CUTADAPT(
        read_pairs_ch,
        sampleseq_file_ch,
        true,  // enable cutadapt
        samples_with_additional_barcodes ?: [:],
        primer_f,
        primer_r,
        false  // skip_polyg_trim
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

    // DADA2 denoising
    dada2_output_ch = WF_DENOISE_DADA2(
        trimmed_sequences_for_dada2,
        dada2_trunc_len_f,
        dada2_trunc_len_r,
        dada2_trim_left_f,
        dada2_trim_left_r,
        1000000,  // n_reads_learn
        0         // n_threads (use available)
    )

    // Filter DADA2 outputs to only include samples present in metadata
    filtering_results_ch = WF_FILTER_BY_METADATA(
        dada2_output_ch.feature_table,
        dada2_output_ch.representative_sequences,
        metadata_ch
    )

    // Run decontamination based on enabled method
    if (use_frequency_decontam) {
        decontam_results_ch = WF_DECONTAM_FREQUENCY(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_ch,
            'DNA concentration'  // default column name
        )
        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
    } else if (use_prevalence_decontam) {
        decontam_results_ch = WF_DECONTAM_PREVALENCE(
            filtering_results_ch.filtered_feature_table,
            filtering_results_ch.filtered_rep_seqs,
            metadata_ch,
            'Decontam Control',     // control column
            'Negative Control'      // control indicator
        )
        filtered_feature_table = decontam_results_ch.filtered_table
        filtered_rep_seqs = decontam_results_ch.filtered_seqs
    } else {
        // Use metadata-filtered results directly
        filtered_feature_table = filtering_results_ch.filtered_feature_table
        filtered_rep_seqs = filtering_results_ch.filtered_rep_seqs
    }

    // Taxonomic Classification with SILVA
    taxonomy_results_ch = CLASSIFY_TAXONOMY(
        filtered_rep_seqs,
        silva_classifier_ch,
        filtered_feature_table,
        metadata_ch
    )

    // Alternative taxonomy classification using DADA2
    dada2_taxonomy_results_ch = WF_DADA2_TAXONOMY(
        filtered_rep_seqs,
        filtered_feature_table
    )

    // Phylogenetic Analysis
    phylogeny_results = ALIGN_AND_TREE_QIIME2(filtered_rep_seqs)

    emit:
    // Core outputs for downstream analysis
    feature_table = filtered_feature_table
    rep_seqs = filtered_rep_seqs
    taxonomy = taxonomy_results_ch.classification
    taxonomy_summary = taxonomy_results_ch.classification_summary
    phylogeny_rooted = phylogeny_results.rooted_tree_qza
    phylogeny_unrooted = phylogeny_results.unrooted_tree_qza
    
    // QC outputs
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    fastqc_results_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.fastqc_results
    multiqc_report_after_cutadapt = WF_QUALITY_CONTROL_AFTER_CUTADAPT.out.multiqc_report
    
    // Denoising outputs
    dada2_denoising_stats = dada2_output_ch.denoising_stats
    dada2_ft_summary = dada2_output_ch.feature_table_summary
    dada2_rep_seqs_summary = dada2_output_ch.rep_seqs_summary
    
    // Alternative taxonomy
    dada2_taxonomy = dada2_taxonomy_results_ch.qiime2_taxonomy
    
    // Classifier
    classifier = silva_classifier_ch
}
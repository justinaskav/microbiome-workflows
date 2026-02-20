// Deblur Denoising Module
// Performs Deblur denoising on paired-end sequences via read joining

// Reuse manifest creation and paired-end import from the DADA2 module
include { CREATE_DADA2_MANIFEST as CREATE_DEBLUR_MANIFEST } from '../dada2/denoising_dada2_wf.nf'
include { IMPORT_ALL_FASTQ_DADA2 as IMPORT_ALL_FASTQ_DEBLUR } from '../dada2/denoising_dada2_wf.nf'

// Reuse summary/tabulation processes with aliases to avoid name collisions
include { SUMMARIZE_FEATURE_TABLE_NO_METADATA as SUMMARIZE_FEATURE_TABLE_DEBLUR } from '../visualization/summary_visualization_wf.nf'
include { TABULATE_SEQS as TABULATE_SEQS_DEBLUR } from '../visualization/summary_visualization_wf.nf'

// Process to join paired-end reads using vsearch
process JOIN_PAIRED_ENDS_DEBLUR {
    tag "join_paired_ends_deblur"
    conda "${params.qiime2_conda_env}"

    input:
    path demux_qza

    output:
    path("joined-seqs.qza"), emit: joined_seqs_qza

    script:
    """
    qiime vsearch merge-pairs \\
        --i-demultiplexed-seqs "${demux_qza}" \\
        --o-merged-sequences joined-seqs.qza \\
        --o-unmerged-sequences unjoined-seqs.qza \\
        --verbose
    """
}

// Process to quality-filter joined sequences before Deblur
process QUALITY_FILTER_DEBLUR {
    tag "quality_filter_deblur"
    conda "${params.qiime2_conda_env}"

    input:
    path joined_seqs_qza

    output:
    path("filtered-seqs.qza"),  emit: filtered_seqs_qza
    path("filter-stats.qza"),   emit: filter_stats_qza

    script:
    """
    qiime quality-filter q-score \\
        --i-demux "${joined_seqs_qza}" \\
        --o-filtered-sequences filtered-seqs.qza \\
        --o-filter-stats filter-stats.qza \\
        --verbose
    """
}

// Process to run Deblur denoise-16S on joined + quality-filtered reads
process RUN_DEBLUR_16S {
    tag "deblur_denoise_16s"
    conda "${params.qiime2_conda_env}"
    label 'process_high'

    input:
    path filtered_seqs_qza
    val trim_length
    val mean_error
    val indel_prob
    val min_reads
    val min_size
    val jobs_to_start

    output:
    path("table.qza"),    emit: feature_table_qza
    path("rep-seqs.qza"), emit: rep_seqs_qza
    path("stats.qza"),    emit: denoising_stats_qza

    script:
    """
    qiime deblur denoise-16S \\
        --i-demultiplexed-seqs "${filtered_seqs_qza}" \\
        --p-trim-length ${trim_length} \\
        --p-mean-error ${mean_error} \\
        --p-indel-prob ${indel_prob} \\
        --p-min-reads ${min_reads} \\
        --p-min-size ${min_size} \\
        --p-jobs-to-start ${jobs_to_start} \\
        --p-sample-stats \\
        --o-table table.qza \\
        --o-representative-sequences rep-seqs.qza \\
        --o-stats stats.qza \\
        --verbose
    """
}

// Process to run Deblur denoise-other with a user-supplied positive-filter reference
process RUN_DEBLUR_OTHER {
    tag "deblur_denoise_other"
    conda "${params.qiime2_conda_env}"
    label 'process_high'

    input:
    path filtered_seqs_qza
    path reference_seqs_qza   // FeatureData[Sequence] — e.g. SILVA V3-V4 extracted reads
    val trim_length
    val mean_error
    val indel_prob
    val min_reads
    val min_size
    val jobs_to_start

    output:
    path("table.qza"),    emit: feature_table_qza
    path("rep-seqs.qza"), emit: rep_seqs_qza
    path("stats.qza"),    emit: denoising_stats_qza

    script:
    """
    qiime deblur denoise-other \\
        --i-demultiplexed-seqs "${filtered_seqs_qza}" \\
        --i-reference-seqs     "${reference_seqs_qza}" \\
        --p-trim-length ${trim_length} \\
        --p-mean-error ${mean_error} \\
        --p-indel-prob ${indel_prob} \\
        --p-min-reads ${min_reads} \\
        --p-min-size ${min_size} \\
        --p-jobs-to-start ${jobs_to_start} \\
        --p-sample-stats \\
        --o-table table.qza \\
        --o-representative-sequences rep-seqs.qza \\
        --o-stats stats.qza \\
        --verbose
    """
}

// Process to visualize per-sample Deblur statistics
process VISUALIZE_DEBLUR_STATS {
    tag "visualize_deblur_stats"
    conda "${params.qiime2_conda_env}"

    input:
    path stats_qza

    output:
    path("denoising_stats_summary.qzv"), emit: denoising_stats_summary_qzv

    script:
    """
    qiime deblur visualize-stats \\
        --i-deblur-stats "${stats_qza}" \\
        --o-visualization denoising_stats_summary.qzv \\
        --verbose
    """
}

// Process to summarize joined sequences (for QC inspection)
process SUMMARIZE_DEMUX_DEBLUR {
    tag "summarize_demux_deblur"
    conda "${params.qiime2_conda_env}"

    input:
    path joined_seqs_qza

    output:
    path("demux_summary.qzv"), emit: demux_summary_qzv

    script:
    """
    qiime demux summarize \\
        --i-data "${joined_seqs_qza}" \\
        --o-visualization demux_summary.qzv \\
        --verbose
    """
}

workflow WF_DENOISE_DEBLUR {
    take:
    read_pairs_ch           // Channel of [meta, [r1, r2]]
    trim_length             // Required: length to trim joined sequences to (e.g. 390 for V3-V4)
    mean_error              // Default: 0.005
    indel_prob              // Default: 0.01
    min_reads               // Default: 10
    min_size                // Default: 2
    jobs_to_start           // Number of parallel jobs
    sample_prefix_to_remove // Optional prefix to remove from sample IDs
    mode                    // '16s' (default) or 'other' (user-supplied reference)
    reference_seqs          // Channel of path to FeatureData[Sequence] QZA (for mode='other')

    main:
    // 1. Build manifest and import as paired-end QZA
    CREATE_DEBLUR_MANIFEST(
        read_pairs_ch.collect(flat: false),
        sample_prefix_to_remove ?: ''
    )
    manifest_ch = CREATE_DEBLUR_MANIFEST.out.manifest_file

    IMPORT_ALL_FASTQ_DEBLUR(manifest_ch)
    demux_qza_ch = IMPORT_ALL_FASTQ_DEBLUR.out.demux_qza

    // 2. Join paired-end reads
    JOIN_PAIRED_ENDS_DEBLUR(demux_qza_ch)
    joined_seqs_ch = JOIN_PAIRED_ENDS_DEBLUR.out.joined_seqs_qza

    // 3. Quality filter joined reads
    QUALITY_FILTER_DEBLUR(joined_seqs_ch)
    filtered_seqs_ch = QUALITY_FILTER_DEBLUR.out.filtered_seqs_qza

    // 4. Deblur denoising — mode-dependent
    if (mode == 'other') {
        RUN_DEBLUR_OTHER(
            filtered_seqs_ch,
            reference_seqs,
            trim_length,
            mean_error,
            indel_prob,
            min_reads,
            min_size,
            jobs_to_start,
        )
        deblur_results_ch = RUN_DEBLUR_OTHER.out
    } else {
        RUN_DEBLUR_16S(
            filtered_seqs_ch,
            trim_length,
            mean_error,
            indel_prob,
            min_reads,
            min_size,
            jobs_to_start,
        )
        deblur_results_ch = RUN_DEBLUR_16S.out
    }

    // 5. Summaries and visualizations
    SUMMARIZE_DEMUX_DEBLUR(joined_seqs_ch)
    SUMMARIZE_FEATURE_TABLE_DEBLUR(deblur_results_ch.feature_table_qza)
    TABULATE_SEQS_DEBLUR(deblur_results_ch.rep_seqs_qza)
    VISUALIZE_DEBLUR_STATS(deblur_results_ch.denoising_stats_qza)

    emit:
    feature_table            = deblur_results_ch.feature_table_qza
    representative_sequences = deblur_results_ch.rep_seqs_qza
    denoising_stats          = deblur_results_ch.denoising_stats_qza
    demux_summary            = SUMMARIZE_DEMUX_DEBLUR.out.demux_summary_qzv
    feature_table_summary    = SUMMARIZE_FEATURE_TABLE_DEBLUR.out.feature_table_summary_qzv
    rep_seqs_summary         = TABULATE_SEQS_DEBLUR.out.rep_seqs_summary_qzv
    denoising_stats_summary  = VISUALIZE_DEBLUR_STATS.out.denoising_stats_summary_qzv
    imported_demux_qza       = demux_qza_ch
}

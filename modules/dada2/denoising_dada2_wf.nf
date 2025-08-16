// Process to create the manifest file for all samples
process CREATE_DADA2_MANIFEST {
    tag "create_dada2_manifest"

    input:
    val read_pair_list // List of [meta, [r1_path, r2_path]]

    output:
    path "manifest_all_samples.tsv", emit: manifest_file

    script:
    // read_pair_list is a List of [SimpleMap(id:...), [Path, Path]]
    // Each Path object's toString() method gives its path.
    def manifest_lines = ["sample-id\tforward-absolute-filepath\treverse-absolute-filepath"]
    read_pair_list.each { item ->
        def sample_id = item[0].split("_")[0] // Extract sample ID from meta.id
        // These paths are from the fromFilePairs channel.
        // Making them explicitly absolute for the manifest is safer for QIIME2.
        def fwd_path = new File(item[1][0].toString()).getAbsolutePath()
        def rev_path = new File(item[1][1].toString()).getAbsolutePath()
        manifest_lines.add("${sample_id}\t${fwd_path}\t${rev_path}")
    }
    """
    printf '%s\\n' '${manifest_lines.join('\n')}' > manifest_all_samples.tsv
    """
}

// Process to import ALL paired-end FASTQ files from a channel into a single QIIME2 QZA artifact
process IMPORT_ALL_FASTQ_DADA2 {
    tag "import_all_samples_for_dada2"
    conda "${params.qiime2_conda_env}"

    input:
    path manifest_file // A single manifest file for all samples

    output:
    path("all-samples-demux-paired-end.qza"), emit: demux_qza

    script:
    """
    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-path ${manifest_file} \\
        --output-path "all-samples-demux-paired-end.qza" \\
        --input-format PairedEndFastqManifestPhred33V2 \\
    """
}

// Process to run DADA2 on paired-end data (all samples in one QZA)
process RUN_DADA2_PAIRED_END_QIIME2 {
    tag "dada2_denoise_paired"
    conda "${params.qiime2_conda_env}"
    label 'process_high'

    input:
    path demux_qza // The demux.qza from IMPORT_ALL_FASTQ_DADA2
    val trunc_len_f
    val trunc_len_r
    val trim_left_f
    val trim_left_r
    val n_reads_learn
    val dada2_n_threads // Explicitly pass desired threads for DADA2 tool

    output:
    path("table.qza"), emit: feature_table_qza
    path("rep-seqs.qza"), emit: rep_seqs_qza
    path("denoising-stats.qza"), emit: denoising_stats_qza

    script:
    // Use dada2_n_threads if provided, else default to task.cpus, or 1 if task.cpus is 0
    def threads_to_use = dada2_n_threads ?: (task.cpus > 0 ? task.cpus : 1)

    """
    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs "${demux_qza}" \\
        --p-trunc-len-f ${trunc_len_f} \\
        --p-trunc-len-r ${trunc_len_r} \\
        --p-trim-left-f ${trim_left_f} \\
        --p-trim-left-r ${trim_left_r} \\
        --p-n-reads-learn ${n_reads_learn} \\
        --p-n-threads ${threads_to_use} \\
        --o-table table.qza \\
        --o-representative-sequences rep-seqs.qza \\
        --o-denoising-stats denoising-stats.qza \\
        --verbose
    """
}

// Process to summarize demultiplexed sequences
process SUMMARIZE_DEMUX_DADA2 {
    tag "summarize_demux_dada2"
    conda "${params.qiime2_conda_env}"

    input:
    path demux_qza

    output:
    path("demux_summary.qzv"), emit: demux_summary_qzv

    script:
    """
    qiime demux summarize \\
        --i-data "${demux_qza}" \\
        --o-visualization demux_summary.qzv \\
        --verbose
    """
}

// Process to tabulate denoising stats (DADA2-specific)
process TABULATE_DENOISING_STATS_DADA2 {
    tag "tabulate_denoising_stats_dada2"
    conda "${params.qiime2_conda_env}"

    input:
    path denoising_stats_qza

    output:
    path("denoising_stats_summary.qzv"), emit: denoising_stats_summary_qzv

    script:
    """
    qiime metadata tabulate \\
        --m-input-file "${denoising_stats_qza}" \\
        --o-visualization denoising_stats_summary.qzv \\
        --verbose
    """
}

// Include reusable summary and visualization processes
include { SUMMARIZE_FEATURE_TABLE_NO_METADATA; TABULATE_SEQS } from '../visualization/summary_visualization_wf.nf'

workflow WF_DENOISE_DADA2 {
    take:
        read_pairs_ch // Channel of [meta, [r1, r2]]
        dada2_trunc_len_f
        dada2_trunc_len_r
        dada2_trim_left_f
        dada2_trim_left_r
        dada2_n_reads_learn
        dada2_n_threads

    main:
        // 1. Prepare manifest from read_pairs_ch
        CREATE_DADA2_MANIFEST (
            read_pairs_ch.collect(flat: false)
        )
        dada2_manifest_ch = CREATE_DADA2_MANIFEST.out.manifest_file

        // 2. Import all reads into a single QZA
        IMPORT_ALL_FASTQ_DADA2 ( dada2_manifest_ch )
        demux_qza_ch = IMPORT_ALL_FASTQ_DADA2.out.demux_qza

        // 3. Run DADA2 denoise-paired
        RUN_DADA2_PAIRED_END_QIIME2 (
            demux_qza_ch,
            dada2_trunc_len_f,
            dada2_trunc_len_r,
            dada2_trim_left_f,
            dada2_trim_left_r,
            dada2_n_reads_learn,
            dada2_n_threads
        )
        dada2_results_ch = RUN_DADA2_PAIRED_END_QIIME2.out

        // 4. Generate Summaries/Visualizations
        SUMMARIZE_DEMUX_DADA2 ( demux_qza_ch )
        SUMMARIZE_FEATURE_TABLE_NO_METADATA ( dada2_results_ch.feature_table_qza )
        TABULATE_SEQS ( dada2_results_ch.rep_seqs_qza )
        TABULATE_DENOISING_STATS_DADA2 ( dada2_results_ch.denoising_stats_qza )

    emit:
        feature_table = dada2_results_ch.feature_table_qza
        representative_sequences = dada2_results_ch.rep_seqs_qza
        denoising_stats = dada2_results_ch.denoising_stats_qza
        demux_summary = SUMMARIZE_DEMUX_DADA2.out.demux_summary_qzv
        feature_table_summary = SUMMARIZE_FEATURE_TABLE_NO_METADATA.out.feature_table_summary_qzv
        rep_seqs_summary = TABULATE_SEQS.out.rep_seqs_summary_qzv
        denoising_stats_summary = TABULATE_DENOISING_STATS_DADA2.out.denoising_stats_summary_qzv
        imported_demux_qza = demux_qza_ch
}


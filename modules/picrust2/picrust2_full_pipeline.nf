// PICRUSt2 Full Pipeline Process
// Runs the complete PICRUSt2 pipeline for functional prediction

nextflow.enable.dsl = 2

process PICRUST2_FULL_PIPELINE {
    tag "picrust2_full_pipeline"
    conda "${params.qiime2_conda_env}"
    label 'process_high'

    input:
    path feature_table_qza
    path rep_seqs_qza
    val picrust2_threads

    output:
    path "picrust2_output/ko_metagenome.qza", emit: ko_metagenome
    path "picrust2_output/ec_metagenome.qza", emit: ec_metagenome
    path "picrust2_output/pathway_abundance.qza", emit: pathway_abundance

    script:
    // Use picrust2_threads if provided, else default to task.cpus, or 1 if task.cpus is 0
    def threads_to_use = picrust2_threads ?: (task.cpus > 0 ? task.cpus : 1)

    """
    # Create output directory
    mkdir -p picrust2_output

    # Run PICRUSt2 full pipeline
    qiime picrust2 full-pipeline \\
        --i-table ${feature_table_qza} \\
        --i-seq ${rep_seqs_qza} \\
        --p-threads ${threads_to_use} \\
        --p-hsp-method ${params.picrust2?.hsp_method ?: 'mp'} \\
        --p-max-nsti ${params.picrust2?.max_nsti ?: 2} \\
        --o-ko-metagenome picrust2_output/ko_metagenome.qza \\
        --o-ec-metagenome picrust2_output/ec_metagenome.qza \\
        --o-pathway-abundance picrust2_output/pathway_abundance.qza \\
        --verbose
    """
}

// Workflow wrapper for PICRUSt2 full pipeline
workflow WF_PICRUST2_FULL_PIPELINE {
    take:
    feature_table_ch    // path to feature table QZA
    rep_seqs_ch         // path to representative sequences QZA
    picrust2_threads    // number of threads to use

    main:
    picrust2_results = PICRUST2_FULL_PIPELINE(
        feature_table_ch,
        rep_seqs_ch,
        picrust2_threads
    )

    emit:
    ko_metagenome = picrust2_results.ko_metagenome
    ec_metagenome = picrust2_results.ec_metagenome
    pathway_abundance = picrust2_results.pathway_abundance
}

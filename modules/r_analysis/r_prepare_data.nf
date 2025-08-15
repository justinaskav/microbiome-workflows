// R Data Preparation Module
nextflow.enable.dsl = 2

process PREPARE_TSE_DATA {
    tag "prepare_tse_data"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    label 'process_medium'

    input:
    path feature_table_qza
    path taxonomy_qza
    path phylogeny_qza
    path metadata_file
    path prepare_script

    output:
    path "tse_data.rds", emit: tse_file
    path "sample_data.rds", emit: sample_data
    path "feature_data.rds", emit: feature_data

    script:
    """
    Rscript ${prepare_script} \
        --feature_table ${feature_table_qza} \
        --taxonomy ${taxonomy_qza} \
        --phylogeny ${phylogeny_qza} \
        --metadata ${metadata_file} \
        --output_dir . \
        --output_prefix tse
    """
}

process CONVERT_QIIME2_TO_R {
    tag "convert_qiime2"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    
    input:
    path qiime2_artifact
    val artifact_type  // 'feature_table', 'taxonomy', 'phylogeny'
    path conversion_script
    
    output:
    path "converted_${artifact_type}.rds", emit: converted_data
    
    script:
    """
    Rscript ${conversion_script} \
        --input ${qiime2_artifact} \
        --type ${artifact_type} \
        --output converted_${artifact_type}.rds
    """
}
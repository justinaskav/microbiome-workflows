// MaAsLin3 differential abundance analysis process for R microbiome workflows

process CALC_MAASLIN3 {
    tag "${tse_type}_maaslin3"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    publishDir "${params.outputDir}/maaslin3", mode: 'copy'
    
    input:
    tuple val(tse_type), path(tse_file), val(config), path(maaslin3_script)

    output:
    tuple val(tse_type), path("*_maaslin3_results.rds"), emit: maaslin3_results
    tuple val(tse_type), path("*_maaslin3_plots.pdf"), emit: maaslin3_plots, optional: true
    path "*.log", emit: log

    script:
    def prv_cut = config.prv_cut ?: 0.1
    def lib_cut = config.lib_cut ?: 1000
    def alpha = config.alpha ?: 0.05
    def p_adj_method = config.p_adj_method ?: "BH"
    def filter_incertae_sedis = config.filter_incertae_sedis ?: true
    def normalization = config.normalization ?: "TSS"
    def transform = config.transform ?: "LOG"
    
    """
    Rscript ${maaslin3_script} \\
        --tse ${tse_file} \\
        --tse_type ${tse_type} \\
        --formula "${config.formula}" \\
        --output_dir . \\
        --output_prefix ${tse_type} \\
        --prv_cut ${prv_cut} \\
        --lib_cut ${lib_cut} \\
        --alpha ${alpha} \\
        --p_adj_method ${p_adj_method} \\
        --filter_incertae_sedis ${filter_incertae_sedis} \\
        --normalization ${normalization} \\
        --transform ${transform} \\
        2>&1 | tee ${tse_type}_maaslin3.log
    """
}
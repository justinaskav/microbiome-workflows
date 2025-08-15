
// MaAsLin3 workflow module
process MAASLIN3 {
    conda params.r_conda_env
    tag "MaAsLin3-${taxa_level}"
    
    input:
    tuple path(feature_table), path(rooted_tree), path(taxonomy), path(metadata), val(formula), val(reference_levels), val(taxa_level)
    path script
    
    output:
    path "*"    // All MaAsLin3 output files directly in the output directory
    
    script:
    def tree_arg = rooted_tree.name != 'OPTIONAL_FILE' ? "--tree ${rooted_tree}" : ""
    def ref_levels_arg = reference_levels ? "--reference_levels '${reference_levels}'" : ""
    
    """
    Rscript ${script} \
        --features ${feature_table} \
        ${tree_arg} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --formula '${formula}' \
        ${ref_levels_arg} \
        --taxa_level '${taxa_level}' \
        --seed '${params.seed}'
    """
}

workflow MAASLIN3_WORKFLOW {
    take:
    combined_inputs_ch    // Channel of tuples containing all inputs

    main:
    def maaslin3_script = file("${moduleDir}/../scripts/maaslin3.R")
    MAASLIN3(combined_inputs_ch, maaslin3_script)

    emit:
    maaslin3_results_ch = MAASLIN3.out
}


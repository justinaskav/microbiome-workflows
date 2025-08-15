// Standardized R Microbiome Analysis Workflow
// Extracted from microbiome monorepo for reusability

nextflow.enable.dsl = 2

// Include required modules
include { PREPARE_TSE_DATA } from '../modules/r_analysis/r_prepare_data.nf'
include { CALC_ALPHA_DIVERSITY } from '../modules/r_analysis/r_alpha_diversity.nf'
include { CALC_BETA_DIVERSITY } from '../modules/r_analysis/r_beta_diversity.nf'
include { CALC_ANCOMBC2 } from '../modules/r_analysis/r_ancombc2.nf'
include { CALC_MAASLIN3 } from '../modules/r_analysis/r_maaslin3.nf'

// Main reusable workflow for R microbiome analysis
workflow R_MICROBIOME_ANALYSIS {
    take:
    feature_table            // QIIME2 feature table (.qza)
    taxonomy                // QIIME2 taxonomy (.qza)
    phylogeny               // QIIME2 phylogenetic tree (.qza)
    metadata                // Metadata file (.tsv)
    formula                 // Statistical formula (e.g., "Group + Gut.Site")
    analysis_config         // Configuration map with analysis parameters

    main:
    // Default analysis configuration
    def default_config = [
        alpha_indices: "shannon,simpson_lambda,faith_diversity",
        alpha_names: "Shannon,Simpson,Faith_PD",
        alpha_niter: 100,
        prv_cut: 0.1,
        lib_cut: 1000,
        alpha: 0.05,
        p_adj_method: "BH",
        filter_incertae_sedis: true,
        normalization: "TSS",
        transform: "LOG"
    ]
    
    // Merge with user config
    def config = default_config + (analysis_config ?: [:])

    // Set up script channels (these will be copied from the microbiome-shared repository)
    prepare_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/prepare_data.R")
    alpha_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/calc_alpha_diversity.R") 
    beta_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/calc_beta_diversity.R")
    ancombc2_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/calc_ancombc2.R")
    maaslin3_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/calc_maaslin3.R")

    // Prepare TreeSummarizedExperiment (TSE) data from QIIME2 artifacts
    tse_results = PREPARE_TSE_DATA(
        feature_table,
        taxonomy, 
        phylogeny,
        metadata,
        prepare_script
    )

    // Create channels for different taxonomic levels
    // Use combine to properly pair each level with the TSE file
    tse_data_ch = Channel.of(['all', 'Phylum', 'Class', 'Order', 'Family', 'Genus'])
        .combine(tse_results.tse_file)

    // Alpha diversity analysis - use working combine approach
    alpha_results = CALC_ALPHA_DIVERSITY(
        tse_data_ch.combine(alpha_script),
        config.alpha_indices,
        config.alpha_names,
        config.alpha_niter
    )

    // Beta diversity analysis - use working combine approach
    beta_input_ch = alpha_results.tse_with_alpha
        .map { tse_file ->
            // Extract tse_type from filename pattern tse_TYPE_with_alpha.rds
            def tse_type = tse_file.getBaseName().replaceAll(/^tse_/, '').replaceAll(/_with_alpha$/, '')
            tuple(tse_type, tse_file)
        }
    
    beta_results = CALC_BETA_DIVERSITY(
        beta_input_ch.combine(beta_script),
        formula
    )

    // Prepare data with configuration for differential abundance
    tse_with_config_ch = tse_data_ch.map { level, tse_file ->
        def level_config = config + [
            formula: formula,
            tse_type: level
        ]
        tuple(level, tse_file, level_config)
    }

    // Differential abundance analysis with ANCOM-BC2 - use working combine approach
    ancombc2_results = CALC_ANCOMBC2(
        tse_with_config_ch.combine(ancombc2_script)
    )

    // Differential abundance analysis with MaAsLin3 - use working combine approach
    maaslin3_results = CALC_MAASLIN3(
        tse_with_config_ch.combine(maaslin3_script)
    )

    emit:
    // Core TSE data for further analysis
    tse_data = tse_results.tse_file
    
    // Diversity analysis results
    tse_with_alpha = alpha_results.tse_with_alpha
    alpha_diversity_results = alpha_results.alpha_results
    beta_diversity_results = beta_results.beta_results
    beta_diversity_plots = beta_results.beta_plots
    
    // Differential abundance results
    ancombc2_results = ancombc2_results.ancombc2_results
    maaslin3_results = maaslin3_results.maaslin3_results
    
    // Summary outputs for reporting
    analysis_summary = Channel.value([
        config: config,
        formula: formula,
        timestamp: new Date().format("yyyy-MM-dd HH:mm:ss")
    ])
}

// Workflow for taxonomic composition analysis
workflow TAXONOMIC_COMPOSITION_ANALYSIS {
    take:
    tse_data                // TSE data file
    metadata               // Metadata file
    grouping_vars          // Variables to group by (e.g., ["Group", "Gut.Site"])
    
    main:
    composition_script = Channel.fromPath("${projectDir}/../microbiome-shared/R/calc_composition.R")
    
    // Calculate taxonomic composition
    composition_results = CALC_TAXONOMIC_COMPOSITION(
        tse_data,
        metadata,
        composition_script,
        grouping_vars
    )
    
    emit:
    composition_plots = composition_results.composition_plots
    composition_data = composition_results.composition_data
}

// Process for taxonomic composition calculation
process CALC_TAXONOMIC_COMPOSITION {
    tag "taxonomic_composition"
    conda "${params.r_conda_env ?: 'r-microbiome'}"
    
    input:
    path tse_file
    path metadata_file
    path composition_script
    val grouping_vars
    
    output:
    path "composition_plots/*", emit: composition_plots
    path "composition_data.rds", emit: composition_data
    
    script:
    def groups = grouping_vars.join(",")
    """
    Rscript ${composition_script} \
        --tse ${tse_file} \
        --metadata ${metadata_file} \
        --output_dir . \
        --grouping_vars "${groups}" \
        --plot_levels "Phylum,Class,Genus" \
        --top_n 10
    """
}
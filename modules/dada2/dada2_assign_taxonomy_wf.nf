// DADA2's assignTaxonomy and addSpecies processes for taxonomic classification

// Include visualization workflow
include { VISUALIZE_TAXONOMY } from '../taxonomy/taxonomy_visualization_wf.nf'

process DOWNLOAD_SILVA_REFERENCES {
    tag "download-silva-refs"
    label 'process_low'
    
    output:
    path("silva_nr99_v138.2_toGenus_trainset.fa.gz"), emit: trainset_file
    path("silva_v138.2_assignSpecies.fa.gz"), emit: species_file
    
    script:
    """
    # Download SILVA reference files for taxonomy assignment
    curl -L -o silva_nr99_v138.2_toGenus_trainset.fa.gz "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"
    curl -L -o silva_v138.2_assignSpecies.fa.gz "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1"
    """
}

process DADA2_ASSIGN_TAXONOMY {
    tag "dada2-assign-taxonomy"
    label 'process_medium'

    conda "${params.r_conda_env}"

    input:
    path representative_seqs
    path trainset_file
    val minBoot // Minimum bootstrap confidence for assigning a taxonomic level
    path script
    
    output:
    path("taxonomy.tsv"), emit: taxonomy_tsv
    path("taxonomy.rds"), emit: taxonomy_rds
    
    script:
    """
    Rscript ${script} \
        ${representative_seqs} \
        ${trainset_file} \
        ${minBoot} \
        ${task.cpus} \
        ${params.seed}
    """
}

process DADA2_ADD_SPECIES {
    tag "dada2-add-species"
    label 'process_medium'

    conda "${params.r_conda_env}"

    input:
    path taxonomy_rds
    path representative_seqs
    path species_file
    val allow_multiple
    val verbose
    path script
    
    output:
    path("taxonomy_with_species.tsv"), emit: taxonomy_with_species_tsv
    path("taxonomy_with_species.rds"), emit: taxonomy_with_species_rds
    
    script:
    // Pass the actual integer value for allow_multiple
    """
    Rscript ${script} \
        ${taxonomy_rds} \
        ${representative_seqs} \
        ${species_file} \
        ${allow_multiple} \
        ${verbose} \
        ${params.seed}
    """
}

process FORMAT_DADA2_TAXONOMY {
    tag "format-dada2-taxonomy"
    label 'process_low'

    conda "${params.r_conda_env}"

    input:
    path taxonomy_tsv      // Original taxonomy TSV from DADA2
    path format_script     // R script to format taxonomy
    
    output:
    path("formatted_taxonomy.tsv"), emit: formatted_tsv
    
    script:
    """
    Rscript ${format_script} \
        ${taxonomy_tsv} \
        formatted_taxonomy.tsv
    """
}

process CONVERT_DADA2_TAXONOMY_TO_QIIME2 {
    tag "convert-dada2-to-qiime2"
    label 'process_low'

    conda "${params.qiime2_conda_env}"

    input:
    path taxonomy_tsv
    path representative_seqs
    
    output:
    path("dada2-taxonomy.qza"), emit: taxonomy_qza
    
    script:
    """
    # First, prepare the taxonomy file in the exact format QIIME2 expects
    # We need to have just the header and data, no metadata header line
    awk 'NR==1 || NR>1' ${taxonomy_tsv} > clean_taxonomy.tsv
    
    # Import the taxonomy into QIIME2 format
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format TSVTaxonomyFormat \
      --input-path clean_taxonomy.tsv \
      --output-path dada2-taxonomy.qza
    """
}



workflow WF_DADA2_TAXONOMY {
    take:
    representative_seqs      // Path to ASV representative sequences
    feature_table            // Path to feature table (optional, for visualization)
    
    main:
    // Download Silva reference files for taxonomy assignment
    silva_refs_ch = DOWNLOAD_SILVA_REFERENCES()
    
    // Set default parameters if not provided in params
    def minBoot = params.classification_confidence ?: 80 // Use the classification confidence for bootstrap threshold
    def allowMultiple = params.dada2_add_species_allow_multiple ?: 3
    def verbose = params.dada2_add_species_verbose ?: true
    
    // Get paths to the R scripts
    def assign_taxonomy_script = file("${moduleDir}/../../scripts/dada2_assign_taxonomy.R")
    def add_species_script = file("${moduleDir}/../../scripts/dada2_add_species.R")
    def format_taxonomy_script = file("${moduleDir}/../../scripts/format_dada2_taxonomy.R")

    // Assign taxonomy using DADA2's assignTaxonomy
    dada2_tax_ch = DADA2_ASSIGN_TAXONOMY(
        representative_seqs,
        silva_refs_ch.trainset_file,
        minBoot,
        assign_taxonomy_script
    )
    
    // Add species using DADA2's addSpecies
    dada2_species_ch = DADA2_ADD_SPECIES(
        dada2_tax_ch.taxonomy_rds,
        representative_seqs,
        silva_refs_ch.species_file,
        allowMultiple,
        verbose,
        add_species_script
    )
    
    // Format taxonomy for QIIME2 compatibility
    formatted_tax_ch = FORMAT_DADA2_TAXONOMY(
        dada2_species_ch.taxonomy_with_species_tsv,
        format_taxonomy_script,
    )
    
    // Convert to QIIME2 format using formatted taxonomy
    qiime2_tax_ch = CONVERT_DADA2_TAXONOMY_TO_QIIME2(
        formatted_tax_ch.formatted_tsv,
        representative_seqs
    )
    
    // Visualize taxonomy if feature table is provided
    // taxonomy_barplot_ch = Channel.empty()
    // if (feature_table) {
    //     // Create barplot visualization using the common visualization workflow
    //     taxonomy_barplot_ch = VISUALIZE_TAXONOMY(
    //         qiime2_tax_ch.taxonomy_qza,
    //         feature_table,
    //         "dada2_taxonomy_barplot"
    //     ).barplot
    // }
    
    emit:
    dada2_taxonomy = dada2_tax_ch.taxonomy_tsv
    dada2_species_taxonomy = dada2_species_ch.taxonomy_with_species_tsv
    dada2_formatted_taxonomy = formatted_tax_ch.formatted_tsv
    qiime2_taxonomy = qiime2_tax_ch.taxonomy_qza
    // taxonomy_barplot = taxonomy_barplot_ch
}

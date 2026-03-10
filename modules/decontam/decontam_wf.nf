// Decontamination workflow using QIIME2's q2-quality-control plugin
// Based on experiment.py implementation for both frequency and prevalence methods

// Process to identify contaminants using frequency method (based on DNA concentration)
process DECONTAM_IDENTIFY_FREQUENCY {
    tag "decontam_identify_frequency"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza    // Feature table (QZA)
    path metadata_tsv         // Metadata file
    val concentration_column  // Column name containing DNA concentration values

    output:
    path("decontam_scores_frequency.qza"), emit: decontam_scores
    path("decontam_scores_frequency.qzv"), emit: decontam_scores_viz

    script:
    """
    # Identify contaminants using frequency method
    qiime quality-control decontam-identify \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${metadata_tsv} \\
        --p-method frequency \\
        --p-freq-concentration-column '${concentration_column}' \\
        --o-decontam-scores decontam_scores_frequency.qza \\
        --verbose

    # Create visualization of decontam scores
    qiime quality-control decontam-score-viz \\
        --i-decontam-scores decontam_scores_frequency.qza \\
        --i-table ${feature_table_qza} \\
        --o-visualization decontam_scores_frequency.qzv \\
        --verbose
    """
}

// Process to identify contaminants using prevalence method (based on control samples)
process DECONTAM_IDENTIFY_PREVALENCE {
    tag "decontam_identify_prevalence"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza    // Feature table (QZA)
    path metadata_tsv         // Metadata file
    val control_column        // Column name identifying control samples
    val control_indicator     // Value in control column that indicates control samples

    output:
    path("decontam_scores_prevalence.qza"), emit: decontam_scores
    path("decontam_scores_prevalence.qzv"), emit: decontam_scores_viz

    script:
    """
    # Identify contaminants using prevalence method
    qiime quality-control decontam-identify \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${metadata_tsv} \\
        --p-method prevalence \\
        --p-prev-control-column '${control_column}' \\
        --p-prev-control-indicator '${control_indicator}' \\
        --o-decontam-scores decontam_scores_prevalence.qza \\
        --verbose

    # Create visualization of decontam scores
    qiime quality-control decontam-score-viz \\
        --i-decontam-scores decontam_scores_prevalence.qza \\
        --i-table ${feature_table_qza} \\
        --o-visualization decontam_scores_prevalence.qzv \\
        --verbose
    """
}

process DECONTAM_IDENTIFY_COMBINED {
    tag "decontam_identify_combined"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza    // Feature table (QZA)
    path metadata_tsv         // Metadata file
    val concentration_column  // Column name containing DNA concentration values
    val control_column        // Column name identifying control samples
    val control_indicator     // Value in control column that indicates control samples

    output:
    path("decontam_scores_combined.qza"), emit: decontam_scores
    path("decontam_scores_combined.qzv"), emit: decontam_scores_viz

    script:
    """
    # Identify contaminants using combined frequency and prevalence methods
    qiime quality-control decontam-identify \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${metadata_tsv} \\
        --p-method combined \\
        --p-freq-concentration-column '${concentration_column}' \\
        --p-prev-control-column '${control_column}' \\
        --p-prev-control-indicator '${control_indicator}' \\
        --o-decontam-scores decontam_scores_combined.qza \\
        --verbose

    # Create visualization of decontam scores
    qiime quality-control decontam-score-viz \\
        --i-decontam-scores decontam_scores_combined.qza \\
        --i-table ${feature_table_qza} \\
        --o-visualization decontam_scores_combined.qzv \\
        --verbose
    """
}

// Process to remove identified contaminants from feature table using feature-table filter
process DECONTAM_FILTER_TABLE {
    tag "decontam_filter_table"
    conda "${params.qiime2_conda_env}"

    input:
    path decontam_scores_qza  // Decontam scores (QZA)
    path feature_table_qza    // Original feature table (QZA)
    val threshold             // P-value threshold for filtering (default: 0.1)

    output:
    path("decontam_filtered_table.qza"), emit: filtered_table

    script:
    def filter_threshold = threshold ?: 0.1
    """
    # Filter feature table based on decontam scores
    # Keep features where p > THRESHOLD OR p IS NULL (non-contaminants)
    qiime feature-table filter-features \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${decontam_scores_qza} \\
        --p-where '[p]>${filter_threshold} OR [p] IS NULL' \\
        --o-filtered-table decontam_filtered_table.qza \\
        --verbose
    """
}

// Process to filter representative sequences based on filtered feature table
process DECONTAM_FILTER_SEQS {
    tag "decontam_filter_seqs"
    conda "${params.qiime2_conda_env}"

    input:
    path rep_seqs_qza         // Representative sequences (QZA)
    path filtered_table_qza   // Filtered feature table (QZA)

    output:
    path("decontam_filtered_seqs.qza"), emit: filtered_seqs

    script:
    """
    # Filter representative sequences to match filtered feature table
    qiime feature-table filter-seqs \\
        --i-data ${rep_seqs_qza} \\
        --i-table ${filtered_table_qza} \\
        --o-filtered-data decontam_filtered_seqs.qza \\
        --verbose
    """
}

// Process to extract contaminant features from feature table (inverse of DECONTAM_FILTER_TABLE)
process DECONTAM_EXTRACT_CONTAMINANT_TABLE {
    tag "extract_contaminant_table"
    conda "${params.qiime2_conda_env}"

    input:
    path decontam_scores_qza  // Decontam scores (QZA)
    path feature_table_qza    // Original feature table (QZA)
    val threshold             // P-value threshold for filtering (default: 0.1)

    output:
    path("contaminant_table.qza"), emit: contaminant_table

    script:
    def filter_threshold = threshold ?: 0.1
    """
    # Extract contaminant features (p <= threshold = contaminant)
    qiime feature-table filter-features \\
        --i-table ${feature_table_qza} \\
        --m-metadata-file ${decontam_scores_qza} \\
        --p-where '[p]<=${filter_threshold}' \\
        --o-filtered-table contaminant_table.qza \\
        --verbose
    """
}

// Process to extract contaminant sequences based on contaminant feature table
process DECONTAM_EXTRACT_CONTAMINANT_SEQS {
    tag "extract_contaminant_seqs"
    conda "${params.qiime2_conda_env}"

    input:
    path rep_seqs_qza           // Representative sequences (QZA)
    path contaminant_table_qza  // Contaminant feature table (QZA)

    output:
    path("contaminant_seqs.qza"), emit: contaminant_seqs

    script:
    """
    # Filter representative sequences to match contaminant feature table
    qiime feature-table filter-seqs \\
        --i-data ${rep_seqs_qza} \\
        --i-table ${contaminant_table_qza} \\
        --o-filtered-data contaminant_seqs.qza \\
        --verbose
    """
}

// Workflow for frequency-based decontamination
workflow WF_DECONTAM_FREQUENCY {
    take:
    feature_table_qza         // Feature table from DADA2
    rep_seqs_qza             // Representative sequences from DADA2
    metadata_tsv             // Metadata file
    concentration_column     // Column name for DNA concentration (optional, defaults to params)
    threshold                // P-value threshold for filtering (optional, defaults to params.decontam_threshold)

    main:
    // Use provided column name or default parameter
    def conc_column = concentration_column ?: params.decontam?.freq_concentration_column
    def filter_threshold = threshold ?: params.decontam?.threshold

    // Identify contaminants using frequency method
    decontam_identify_results = DECONTAM_IDENTIFY_FREQUENCY(
        feature_table_qza,
        metadata_tsv,
        conc_column
    )

    // Filter feature table based on decontam scores
    filtered_table_results = DECONTAM_FILTER_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Filter representative sequences based on filtered table
    filtered_seqs_results = DECONTAM_FILTER_SEQS(
        rep_seqs_qza,
        filtered_table_results.filtered_table
    )

    // Extract contaminant features (inverse filter)
    contaminant_table_results = DECONTAM_EXTRACT_CONTAMINANT_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Extract contaminant sequences
    contaminant_seqs_results = DECONTAM_EXTRACT_CONTAMINANT_SEQS(
        rep_seqs_qza,
        contaminant_table_results.contaminant_table
    )

    emit:
    decontam_scores = decontam_identify_results.decontam_scores
    decontam_scores_viz = decontam_identify_results.decontam_scores_viz
    filtered_table = filtered_table_results.filtered_table
    filtered_seqs = filtered_seqs_results.filtered_seqs
    contaminant_table = contaminant_table_results.contaminant_table
    contaminant_seqs = contaminant_seqs_results.contaminant_seqs
}

// Workflow for prevalence-based decontamination
workflow WF_DECONTAM_PREVALENCE {
    take:
    feature_table_qza         // Feature table from DADA2
    rep_seqs_qza             // Representative sequences from DADA2
    metadata_tsv             // Metadata file
    control_column           // Column name for control identification (optional, defaults to params)
    control_indicator        // Value indicating control samples (optional, defaults to params)
    threshold                // P-value threshold for filtering (optional, defaults to params.decontam_threshold)

    main:
    // Use provided parameters or defaults
    def ctrl_column = control_column ?: params.decontam?.prev_control_column
    def ctrl_indicator = control_indicator ?: params.decontam?.prev_control_indicator
    def filter_threshold = threshold ?: params.decontam?.threshold

    // Identify contaminants using prevalence method
    decontam_identify_results = DECONTAM_IDENTIFY_PREVALENCE(
        feature_table_qza,
        metadata_tsv,
        ctrl_column,
        ctrl_indicator
    )

    // Filter feature table based on decontam scores
    filtered_table_results = DECONTAM_FILTER_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Filter representative sequences based on filtered table
    filtered_seqs_results = DECONTAM_FILTER_SEQS(
        rep_seqs_qza,
        filtered_table_results.filtered_table
    )

    // Extract contaminant features (inverse filter)
    contaminant_table_results = DECONTAM_EXTRACT_CONTAMINANT_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Extract contaminant sequences
    contaminant_seqs_results = DECONTAM_EXTRACT_CONTAMINANT_SEQS(
        rep_seqs_qza,
        contaminant_table_results.contaminant_table
    )

    emit:
    decontam_scores = decontam_identify_results.decontam_scores
    decontam_scores_viz = decontam_identify_results.decontam_scores_viz
    filtered_table = filtered_table_results.filtered_table
    filtered_seqs = filtered_seqs_results.filtered_seqs
    contaminant_table = contaminant_table_results.contaminant_table
    contaminant_seqs = contaminant_seqs_results.contaminant_seqs
}

workflow WF_DECONTAM_COMBINED {
    take:
    feature_table_qza         // Feature table from DADA2
    rep_seqs_qza             // Representative sequences from DADA2
    metadata_tsv             // Metadata file
    concentration_column     // Column name for DNA concentration (optional, defaults to params)
    control_column           // Column name for control identification (optional, defaults to params)
    control_indicator        // Value indicating control samples (optional, defaults to params)
    threshold                // P-value threshold for filtering (optional, defaults to params.decontam_threshold)

    main:
    // Use provided parameters or defaults
    def conc_column = concentration_column ?: params.decontam?.freq_concentration_column
    def ctrl_column = control_column ?: params.decontam?.prev_control_column
    def ctrl_indicator = control_indicator ?: params.decontam?.prev_control_indicator
    def filter_threshold = threshold ?: params.decontam?.threshold

    // Identify contaminants using combined method
    decontam_identify_results = DECONTAM_IDENTIFY_COMBINED(
        feature_table_qza,
        metadata_tsv,
        conc_column,
        ctrl_column,
        ctrl_indicator
    )

    // Filter feature table based on decontam scores
    filtered_table_results = DECONTAM_FILTER_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Filter representative sequences based on filtered table
    filtered_seqs_results = DECONTAM_FILTER_SEQS(
        rep_seqs_qza,
        filtered_table_results.filtered_table
    )

    // Extract contaminant features (inverse filter)
    contaminant_table_results = DECONTAM_EXTRACT_CONTAMINANT_TABLE(
        decontam_identify_results.decontam_scores,
        feature_table_qza,
        filter_threshold
    )

    // Extract contaminant sequences
    contaminant_seqs_results = DECONTAM_EXTRACT_CONTAMINANT_SEQS(
        rep_seqs_qza,
        contaminant_table_results.contaminant_table
    )

    emit:
    decontam_scores = decontam_identify_results.decontam_scores
    decontam_scores_viz = decontam_identify_results.decontam_scores_viz
    filtered_table = filtered_table_results.filtered_table
    filtered_seqs = filtered_seqs_results.filtered_seqs
    contaminant_table = contaminant_table_results.contaminant_table
    contaminant_seqs = contaminant_seqs_results.contaminant_seqs
}

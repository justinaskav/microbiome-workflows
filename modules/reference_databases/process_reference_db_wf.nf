// Unified reference database processing workflow
// Handles SILVA and UNITE specific processing steps with parameterized logic

process REVERSE_TRANSCRIBE_RNA {
    tag "$meta.id - RNA to DNA"
    label 'process_medium'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(rna_sequences_qza)

    output:
    tuple val(meta), path("${meta.id}-dna-seqs.qza"), emit: dna_sequences


    script:
    """
    qiime rescript reverse-transcribe \
        --i-rna-sequences "${rna_sequences_qza}" \
        --o-dna-sequences "${meta.id}-dna-seqs.qza" \
        --verbose
    """
}

process CULL_SEQUENCES {
    tag "$meta.id - Cull sequences"
    label 'process_medium'
    cpus 4

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza)

    output:
    tuple val(meta), path("${meta.id}-culled.qza"), emit: culled_sequences


    script:
    """
    qiime rescript cull-seqs \
        --i-sequences "${sequences_qza}" \
        --p-n-jobs ${task.cpus} \
        --o-clean-sequences "${meta.id}-culled.qza" \
        --verbose
    """
}

process FILTER_SEQUENCES_BY_LENGTH {
    tag "$meta.id - Length filter"
    label 'process_medium'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val filter_labels_list
    val filter_min_lengths_list

    output:
    tuple val(meta), path("${meta.id}-len-filt.qza"), emit: filtered_sequences


    script:
    def labels_str = filter_labels_list.join(' ')
    def min_lens_str = filter_min_lengths_list.join(' ')
    """
    qiime rescript filter-seqs-length-by-taxon \
        --i-sequences "${sequences_qza}" \
        --i-taxonomy "${taxonomy_qza}" \
        --p-labels ${labels_str} \
        --p-min-lens ${min_lens_str} \
        --o-filtered-seqs "${meta.id}-len-filt.qza" \
        --o-discarded-seqs "${meta.id}-len-discard.qza" \
        --verbose
    """
}

process EDIT_TAXONOMY {
    tag "$meta.id - Edit taxonomy"
    label 'process_low'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(taxonomy_qza)
    val search_strings
    val replacement_strings

    output:
    tuple val(meta), path("${meta.id}-edited-tax.qza"), emit: edited_taxonomy


    script:
    def search_str = search_strings.join(' ')
    def replace_str = replacement_strings.join(' ')
    """
    qiime rescript edit-taxonomy \
        --i-taxonomy "${taxonomy_qza}" \
        --o-edited-taxonomy "${meta.id}-edited-tax.qza" \
        --p-search-strings '${search_str}' \
        --p-replacement-strings '${replace_str}' \
        --p-use-regex
    """
}

process DEREPLICATE_DATABASE {
    tag "$meta.id - Dereplicate"
    label 'process_high'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val derep_mode

    output:
    tuple val(meta), path("${meta.id}-seqs-derep.qza"), path("${meta.id}-tax-derep.qza"), emit: dereplicated_data

    script:
    """
    qiime rescript dereplicate \
        --i-sequences "${sequences_qza}" \
        --i-taxa "${taxonomy_qza}" \
        --p-mode '${derep_mode}' \
        --p-threads ${task.cpus} \
        --o-dereplicated-sequences "${meta.id}-seqs-derep.qza" \
        --o-dereplicated-taxa "${meta.id}-tax-derep.qza" \
        --verbose
    """
}

process EXTRACT_READS {
    tag "$meta.id - Extract reads"
    label 'process_medium'
    cpus 4

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza)
    val f_primer
    val r_primer

    output:
    tuple val(meta), path("${meta.id}-extracted.qza"), emit: extracted_sequences

    script:
    """
    qiime feature-classifier extract-reads \
        --i-sequences "${sequences_qza}" \
        --p-f-primer '${f_primer}' \
        --p-r-primer '${r_primer}' \
        --p-n-jobs ${task.cpus} \
        --p-read-orientation 'forward' \
        --o-reads "${meta.id}-extracted.qza" \
        --verbose
    """
}

include { DOWNLOAD_REFERENCE_DATABASE } from './download_reference_db_wf.nf'

workflow WF_PROCESS_REF_DB {
    take:
    meta
    db_type              // 'silva' or 'unite'
    db_config            // Map with database-specific parameters
    processing_config    // Map with processing parameters

    main:
    // Download the database
    download_result = DOWNLOAD_REFERENCE_DATABASE(meta, db_type, db_config)
    
    // Initialize processing channel
    current_seqs_ch = download_result.reference_data.map { _meta, seqs, tax -> tuple(_meta, seqs) }
    current_tax_ch = download_result.reference_data.map { _meta, seqs, tax -> tuple(_meta, tax) }

    // SILVA-specific processing
    if (db_type == 'silva') {
        // Reverse transcribe RNA to DNA
        rt_seqs_ch = REVERSE_TRANSCRIBE_RNA(current_seqs_ch).dna_sequences
        current_seqs_ch = rt_seqs_ch

        // Cull sequences
        if (processing_config.cull_sequences ?: true) {
            culled_ch = CULL_SEQUENCES(current_seqs_ch).culled_sequences
            current_seqs_ch = culled_ch
        }

        // Length filtering (SILVA specific)
        if (processing_config.filter_by_length && processing_config.filter_labels && processing_config.filter_min_lengths) {
            filter_input_ch = current_seqs_ch.join(current_tax_ch, by: 0)
            filtered_ch = FILTER_SEQUENCES_BY_LENGTH(
                filter_input_ch,
                processing_config.filter_labels,
                processing_config.filter_min_lengths
            ).filtered_sequences
            current_seqs_ch = filtered_ch
        }
    }

    // UNITE-specific processing  
    if (db_type == 'unite') {
        // Edit taxonomy to remove SH identifiers
        if (processing_config.remove_sh_identifiers ?: true) {
            tax_edit_ch = EDIT_TAXONOMY(
                current_tax_ch,
                [';sh__.*'],
                ['']
            ).edited_taxonomy
            current_tax_ch = tax_edit_ch
        }
    }

    // Common processing steps
    // Prepare for dereplication
    derep_input_ch = current_seqs_ch.join(current_tax_ch, by: 0)
    
    // Dereplicate
    dereplicated_ch = DEREPLICATE_DATABASE(
        derep_input_ch,
        processing_config.derep_mode ?: 'uniq'
    ).dereplicated_data

    // Extract reads if primers provided (common for both databases)
    if (processing_config.f_primer && processing_config.r_primer) {
        extracted_seqs_ch = EXTRACT_READS(
            dereplicated_ch.map { _meta, seqs, tax -> tuple(_meta, seqs) },
            processing_config.f_primer,
            processing_config.r_primer
        ).extracted_sequences

        // Rejoin with taxonomy
        final_seqs_ch = extracted_seqs_ch
        final_tax_ch = dereplicated_ch.map { _meta, seqs, tax -> tuple(_meta, tax) }
    } else {
        // Use full-length dereplicated sequences
        final_seqs_ch = dereplicated_ch.map { _meta, seqs, tax -> tuple(_meta, seqs) }
        final_tax_ch = dereplicated_ch.map { _meta, seqs, tax -> tuple(_meta, tax) }
    }

    // Combine final sequences and taxonomy
    processed_data_ch = final_seqs_ch.join(final_tax_ch, by: 0)

    emit:
    processed_data = processed_data_ch  // [meta, sequences.qza, taxonomy.qza]
}

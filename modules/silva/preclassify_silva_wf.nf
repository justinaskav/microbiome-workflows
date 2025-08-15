process REVERSE_TRANSCRIBE_SILVA {
    tag "$meta.id"
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

process CULL_SEQS_SILVA {
    tag "$meta.id"
    label 'process_medium'
    cpus 4

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(dna_sequences_qza)

    output:
    // Because it will be joined by meta
    tuple val(meta), path("${meta.id}-culled.qza"), emit: culled_sequences

    script:
    """
    qiime rescript cull-seqs \
        --i-sequences "${dna_sequences_qza}" \
        --p-n-jobs ${task.cpus} \
        --o-clean-sequences "${meta.id}-culled.qza" \
        --verbose
    """
}

process FILTER_SEQS_LENGTH_BY_TAXON_SILVA {
    tag "$meta.id"
    label 'process_medium'

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val filter_labels_list
    val filter_min_lengths_list

    output:
    tuple val(meta), path("${meta.id}-len-filt.qza"), emit: filtered_sequences
    // tuple val(meta), path("${meta.id}-len-discard.qza"), emit: discarded_sequences // Optional: uncomment to get discarded sequences
    
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

process DEREPLICATE_SILVA {
    tag "$meta.id"
    label 'process_high' // Can be memory intensive
    cpus 2

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(sequences_qza), path(taxonomy_qza)
    val derep_mode_str

    output:
    tuple val(meta), path("${meta.id}-seqs-derep.qza"), path("${meta.id}-tax-derep.qza"), emit: dereplicated_data

    script:
    """
    qiime rescript dereplicate \
        --i-sequences "${sequences_qza}" \
        --i-taxa "${taxonomy_qza}" \
        --p-mode '${derep_mode_str}' \
        --p-threads ${task.cpus} \
        --o-dereplicated-sequences "${meta.id}-seqs-derep.qza" \
        --o-dereplicated-taxa "${meta.id}-tax-derep.qza" \
        --verbose
    """
}

process EXTRACT_AMPLICON_REGION {
    tag "$meta.id (Amplicon)"
    label 'process_medium'

    conda "${params.qiime2_conda_env}"
    cpus 4

    input:
    tuple val(meta), path(sequences_qza)
    val f_primer_seq
    val r_primer_seq

    output:
    tuple val(meta), path("${meta.id}-amplicon-reads.qza"), emit: amplicon_reads

    script:
    // Workflow logic ensures this process only runs if primers are provided and non-empty.
    """
    qiime feature-classifier extract-reads \
        --i-sequences "${sequences_qza}" \
        --p-f-primer '${f_primer_seq}' \
        --p-r-primer '${r_primer_seq}' \
        --p-n-jobs ${task.cpus} \
        --p-read-orientation 'forward' \
        --o-reads "${meta.id}-amplicon-reads.qza" \
        --verbose
    """
}

process CREATE_CLASSIFIER_QZA {
    tag "$meta.id (Classifier)"
    label 'process_high' 

    conda "${params.qiime2_conda_env}"

    input:
    tuple val(meta), path(reference_reads_qza), path(reference_taxonomy_qza)

    output:
    path("${meta.id}-classifier.qza"), emit: classifier

    script:
    // NOTE: Ngram range is [8, 8] in DADA2, we'll try to use the same.
    """
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "${reference_reads_qza}" \
        --i-reference-taxonomy "${reference_taxonomy_qza}" \
        --p-feat-ext--ngram-range '[8, 8]' \
        --o-classifier "${meta.id}-classifier.qza" \
        --verbose
    """
}

// Import the DOWNLOAD_SILVA workflow
include { DOWNLOAD_SILVA } from './download_silva_wf.nf'

workflow PRECLASSIFY_SILVA {
    take:
    meta
    silva_version
    silva_target
    filter_labels
    filter_min_lengths
    derep_mode
    f_primer
    r_primer

    main:
    // Get SILVA data and create a channel with RNA sequences and taxonomy

    def include_species_labels = false // Set to true if you want to include species labels

    raw_silva_genus_ch = DOWNLOAD_SILVA(
        meta,
        silva_version,
        silva_target,
        include_species_labels
    )

    // Prepare input for reverse transcription by creating a meta map
    rt_input_ch = raw_silva_genus_ch.map { _meta, rna_seqs, _tax_qza -> 
        tuple(_meta, rna_seqs)
    }
    
    // Convert RNA to DNA sequences
    dna_seqs_ch = REVERSE_TRANSCRIBE_SILVA(rt_input_ch).dna_sequences

    // Channel emitting: [meta, culled_seqs.qza]
    culled_seqs_ch = CULL_SEQS_SILVA ( dna_seqs_ch ).culled_sequences

    // Channel emitting: [meta, tax.qza] - for joining with processed sequences
    original_tax_ch = raw_silva_genus_ch.map { _meta, _rna_seqs_path, tax_qza_path -> tuple(_meta, tax_qza_path) }

    // Prepare input for FILTER_SEQS_LENGTH_BY_TAXON_SILVA: [meta, culled_seqs.qza, tax.qza]
    filter_length_input_ch = culled_seqs_ch.join(original_tax_ch, by: 0) // Join by meta (first element)

    // Channel emitting: [meta, len_filt.qza]
    filtered_seqs_len_ch = FILTER_SEQS_LENGTH_BY_TAXON_SILVA (
        filter_length_input_ch,
        filter_labels,
        filter_min_lengths
    ).filtered_sequences

    // Prepare input for DEREPLICATE_SILVA: [meta, len_filt.qza, tax.qza]
    dereplicate_input_ch = filtered_seqs_len_ch.join(original_tax_ch, by: 0) // Join by meta

    // Channel emitting: [meta, seqs_derep.qza, tax_derep.qza]
    dereplicate_silva_ch = DEREPLICATE_SILVA (
        dereplicate_input_ch,
        derep_mode
    ).dereplicated_data

    if (f_primer && r_primer) { // True if primers are non-empty strings
        // Input for EXTRACT_AMPLICON_REGION: [meta, seqs_derep.qza]
        extract_amplicon_input_ch = dereplicate_silva_ch.map { _meta, seqs_derep, _tax_derep -> tuple(_meta, seqs_derep) }

        // Channel emitting: [meta, amplicon_reads.qza]
        amplicon_reads_ch = EXTRACT_AMPLICON_REGION (
            extract_amplicon_input_ch,
            f_primer,
            r_primer
        ).amplicon_reads

        // Need to join amplicon reads with their corresponding DEREPLICATED taxonomy
        // Channel emitting: [meta, tax_derep.qza]
        dereplicated_tax_for_join_ch = dereplicate_silva_ch.map { _meta, _seqs, tax -> tuple(_meta, tax) }
        
        // Final input for CREATE_CLASSIFIER_QZA: [meta, amplicon_reads.qza, tax_derep.qza]
        classifier_input_ch = amplicon_reads_ch.join(dereplicated_tax_for_join_ch, by: 0) // Join by meta
    } else {
        // Use full-length dereplicated sequences and their corresponding dereplicated taxonomy
        // Channel emitting: [meta, seqs_derep.qza, tax_derep.qza]
        classifier_input_ch = dereplicate_silva_ch 
    }

    classifier = CREATE_CLASSIFIER_QZA ( classifier_input_ch ).classifier

    emit:
        classifier
}

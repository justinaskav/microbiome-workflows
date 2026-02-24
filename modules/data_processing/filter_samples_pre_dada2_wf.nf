
// Pre-DADA2 sample filtering: keep only read pairs whose sample ID
// appears in the metadata file.  This avoids wasting compute on samples
// that are not part of the study design.

nextflow.enable.dsl = 2

// Extract sample IDs from a QIIME2-format metadata TSV.
// Skips comment / directive lines (#q2:types, etc.) and outputs one ID per line.
process EXTRACT_METADATA_SAMPLE_IDS {
    tag "extract_metadata_sample_ids"

    input:
    path metadata_file

    output:
    path("sample_ids.txt"), emit: sample_ids

    script:
    """
    awk -F'\\t' '
        NR == 1 { next }        # skip header
        /^#/   { next }        # skip #q2:types and any other directives
        NF > 0 { print \$1 }   # first column = SampleID
    ' ${metadata_file} > sample_ids.txt
    """
}

// Pass-through process that emits reads only when the sample is in the
// metadata list.  The sample ID list is passed as a val so that no file
// staging is needed in the when: block.
process FILTER_READ_PAIR {
    tag "filter_pre_dada2_${meta}"
    cache 'lenient'

    input:
    tuple val(meta), path(reads)
    val sample_ids_list

    output:
    tuple val(meta), path(reads), emit: kept_reads

    when:
    sample_ids_list.contains(meta)

    script:
    """
    true
    """
}

// Workflow: filter a channel of [meta, [r1, r2]] read pairs so that only
// samples listed in the metadata file are retained.
workflow WF_FILTER_SAMPLES_PRE_DADA2 {
    take:
    read_pairs_ch   // channel of [meta, [r1, r2]]
    metadata_file   // path to QIIME2-format metadata TSV

    main:

    // 1. Extract sample IDs from metadata into a file, then read into a list
    EXTRACT_METADATA_SAMPLE_IDS(metadata_file)

    sample_ids_list = EXTRACT_METADATA_SAMPLE_IDS.out.sample_ids
        .splitText()
        .map { it.trim() }
        .filter { it }
        .collect()

    // 2. Filter read pairs using the when: directive in the process
    FILTER_READ_PAIR(
        read_pairs_ch,
        sample_ids_list,
    )

    filtered_reads = FILTER_READ_PAIR.out.kept_reads

    emit:
    filtered_reads
}

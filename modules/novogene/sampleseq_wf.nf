// Novogene-specific SampleSeq_info.xls parsing workflow

process PARSE_SAMPLESEQ_INFO {
    tag "Parse SampleSeq_info.xls"

    input:
    path sampleseq_file

    output:
    path "sample_barcode_primer_info.tsv", emit: barcode_primer_info

    script:
    """
    python3 ${moduleDir}/../../scripts/parse_sampleseq_info.py \
    ${sampleseq_file} \
    -o sample_barcode_primer_info.tsv
    """
}

// Workflow for parsing Novogene SampleSeq info
workflow WF_PARSE_NOVOGENE_SAMPLESEQ {
    take:
    sampleseq_file_ch

    main:
    barcode_primer_info_ch = PARSE_SAMPLESEQ_INFO(
        sampleseq_file_ch,
    ).barcode_primer_info

    emit:
    barcode_primer_info = barcode_primer_info_ch
}

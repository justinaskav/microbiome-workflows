// Novogene-specific SampleSeq_info.xls parsing workflow

process PARSE_SAMPLESEQ_INFO {
    tag "Parse SampleSeq_info.xls"

    input:
    path sampleseq_file
    val samples_with_additional_barcodes

    output:
    path "sample_barcode_primer_info.tsv", emit: barcode_primer_info

    script:
    def additional_barcodes_cmd = ""
    if (samples_with_additional_barcodes) {
        additional_barcodes_cmd = samples_with_additional_barcodes.collect { sample, barcodes ->
            "${sample}\t${barcodes[0]}\t${barcodes[1]}"
        }.join('\n')
    }
    if (additional_barcodes_cmd) {
        """
        python3 ${moduleDir}/../../scripts/parse_sampleseq_info.py \
            ${sampleseq_file} \
            --additional-barcodes "${additional_barcodes_cmd}" \
            -o sample_barcode_primer_info.tsv
        """
    } else {
        """
        python3 ${moduleDir}/../../scripts/parse_sampleseq_info.py \
            ${sampleseq_file} \
            -o sample_barcode_primer_info.tsv
        """
    }
}

// Workflow for parsing Novogene SampleSeq info
workflow WF_PARSE_NOVOGENE_SAMPLESEQ {
    take:
    sampleseq_file_ch
    samples_with_additional_barcodes

    main:
    barcode_primer_info_ch = PARSE_SAMPLESEQ_INFO(
        sampleseq_file_ch,
        samples_with_additional_barcodes
    ).barcode_primer_info

    emit:
    barcode_primer_info = barcode_primer_info_ch
}
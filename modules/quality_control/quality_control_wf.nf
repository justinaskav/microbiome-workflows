// Process to run FastQC on read pairs
process RUN_FASTQC {
    cache 'lenient'
    conda 'bioconda::fastqc'
    cpus params.qc?.fastqc_cpus ?: 2

    input:
    tuple val(meta), path(reads) // Expects e.g. [ id:"sample1", [read1.fq.gz, read2.fq.gz] ]

    output:
    tuple val(meta), path("*.{html,zip}"), emit: fastqc_results // Emits [ meta, [sample_R1_fastqc.html, ..., sample_R2_fastqc.zip] ]

    script:
    def threads = task.cpus > 0 ? task.cpus : 1
    """
    fastqc -t ${threads} -o . --noextract ${reads[0]} ${reads[1]}
    """
}

// Process to run MultiQC to aggregate FastQC reports
process RUN_MULTIQC {
    conda 'bioconda::multiqc python=3.12'

    input:
    path fastqc_dirs // Expects a list of paths from FastQC output (specifically the directories or relevant files for multiqc)  

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data_dir // Optional: MultiQC data directory

    script:
    """
    multiqc . -o .
    """
}


// Workflow within the quality_control_module.nf
workflow WF_QUALITY_CONTROL {
    take:
        read_pairs_ch // input channel: e.g., [ [id:"sample1", [r1,r2]], [id:"sample2", [r1,r2]] ]

    main:
        // Run FastQC for each sample
        fastqc_ch = RUN_FASTQC(read_pairs_ch)

        // Collect all FastQC output directories/files for MultiQC
        multiqc_input_ch = fastqc_ch.fastqc_results.map { _meta, files -> files }.collect()

        // Run MultiQC
        multiqc_ch = RUN_MULTIQC(multiqc_input_ch)

    emit:
        multiqc_report = multiqc_ch.report
        multiqc_data_dir = multiqc_ch.data_dir
        fastqc_results = fastqc_ch.fastqc_results
}

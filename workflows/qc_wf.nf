// Reusable Quality Control Workflow
// All parameters should be set in project-specific nextflow.config files

nextflow.enable.dsl = 2

include { WF_QUALITY_CONTROL } from '../modules/quality_control/quality_control_wf.nf'

// Named workflow for quality control only
workflow QC_ONLY {
    take:
    reads_pattern       // File pattern for paired-end reads

    main:
    read_pairs_ch = Channel.fromFilePairs(reads_pattern)
    WF_QUALITY_CONTROL(read_pairs_ch)

    emit:
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    multiqc_data_dir = WF_QUALITY_CONTROL.out.multiqc_data_dir
}
// Quality Control Subworkflow
// Runs FastQC and MultiQC on read pairs

nextflow.enable.dsl = 2

include { WF_QUALITY_CONTROL } from '../modules/quality_control/quality_control_wf.nf'

workflow QC {
    take:
    read_pairs_ch  // Channel of [sample_id, [r1, r2]]

    main:
    WF_QUALITY_CONTROL(read_pairs_ch)

    emit:
    fastqc_results = WF_QUALITY_CONTROL.out.fastqc_results
    multiqc_report = WF_QUALITY_CONTROL.out.multiqc_report
    multiqc_data_dir = WF_QUALITY_CONTROL.out.multiqc_data_dir
}

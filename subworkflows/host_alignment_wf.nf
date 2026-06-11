// Host Alignment Subworkflow
// Aligns trimmed reads against a prebuilt host bowtie2 index (e.g. GRCh38 + chrM)
// to (a) produce a per-sample MultiQC report of host mapping rates and
// (b) optionally drop host-mapped read pairs before DADA2.
//
// The index is built once at the top level via WF_PREPARE_HOST_INDEX and
// threaded in here, so per-batch invocations only run alignment + MultiQC.
//
// Placement: must run AFTER cutadapt — bowtie2's default end-to-end mode
// requires the full read to align, so primers should already be trimmed.
//
// reads_out:
//   - remove_aligned = true  -> host-unmapped concordant pairs
//   - remove_aligned = false -> input read pairs unchanged (report-only)
//   - enabled = false        -> input read pairs unchanged (no alignment)

nextflow.enable.dsl = 2

include { WF_BOWTIE2_HOST_ALIGNMENT } from '../modules/host_alignment/bowtie2_host_wf.nf'

workflow HOST_ALIGNMENT {
    take:
    read_pairs_ch    // [sample_id, [r1, r2]]
    host_options     // Map: enabled, remove_aligned, sensitivity, cpus
    host_index_ch    // Prebuilt bowtie2 index channel from WF_PREPARE_HOST_INDEX

    main:
    def opts = host_options ?: [:]
    def enabled = opts.enabled != null
        ? opts.enabled
        : (params.host_alignment?.enabled ?: false)
    def remove_aligned = opts.remove_aligned != null
        ? opts.remove_aligned
        : (params.host_alignment?.remove_aligned != null
            ? params.host_alignment.remove_aligned
            : true)

    if (enabled) {
        WF_BOWTIE2_HOST_ALIGNMENT(read_pairs_ch, host_index_ch)
        multiqc_report_ch = WF_BOWTIE2_HOST_ALIGNMENT.out.multiqc_report
        logs_ch           = WF_BOWTIE2_HOST_ALIGNMENT.out.logs
        reads_out_ch      = remove_aligned
            ? WF_BOWTIE2_HOST_ALIGNMENT.out.unmapped_reads
            : read_pairs_ch
    } else {
        multiqc_report_ch = Channel.empty()
        logs_ch           = Channel.empty()
        reads_out_ch      = read_pairs_ch
    }

    emit:
    reads_out      = reads_out_ch
    multiqc_report = multiqc_report_ch
    logs           = logs_ch
}

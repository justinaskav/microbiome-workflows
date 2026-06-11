// Host alignment QC + removal via bowtie2
// Aligns trimmed reads against a host reference (e.g. GRCh38 + chrM), emits
// per-sample bowtie2 logs plus an aggregated MultiQC report, and writes the
// host-unmapped concordant pairs as gzipped FASTQs ready for DADA2.
//
// host_fasta input may be either a local path or a URL (http/https/ftp).
// The subworkflow WF_PREPARE_HOST_INDEX routes URLs through DOWNLOAD_HOST_FASTA.
//
// Designed to run AFTER cutadapt: bowtie2 default end-to-end mode requires
// the full read to align, so primers should already be removed.

nextflow.enable.dsl = 2

process DOWNLOAD_HOST_FASTA {
    tag "Download host reference"
    label 'process_low'

    conda 'conda-forge::curl'

    input:
    val url

    output:
    path "host_reference.fa.gz", emit: fasta

    script:
    """
    curl -L --fail -o host_reference.fa.gz "${url}"
    """
}

process BUILD_BOWTIE2_INDEX {
    tag "Build bowtie2 host index"
    label 'process_medium'
    cache 'lenient'

    conda 'bioconda::bowtie2 bioconda::samtools'

    input:
    path host_fasta

    output:
    path "host_index*.bt2*", emit: index_files

    script:
    def threads = task.cpus > 0 ? task.cpus : 4
    """
    if [[ "${host_fasta}" == *.gz ]]; then
        gunzip -c ${host_fasta} > host_reference.fa
        bowtie2-build --threads ${threads} host_reference.fa host_index
        rm host_reference.fa
    else
        bowtie2-build --threads ${threads} ${host_fasta} host_index
    fi
    """
}

process ALIGN_HOST_BOWTIE2 {
    tag "Host align ${sample_id}"
    cache 'lenient'
    cpus params.host_alignment?.cpus ?: 4

    conda 'bioconda::bowtie2 bioconda::samtools'

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    path "${sample_id}.bowtie2.log", emit: log
    tuple val(sample_id),
          path("${sample_id}_unmapped_R1.fastq.gz"),
          path("${sample_id}_unmapped_R2.fastq.gz"), emit: unmapped

    script:
    def sensitivity = params.host_alignment?.sensitivity ?: 'very-sensitive'
    def threads = task.cpus > 0 ? task.cpus : 4
    """
    bowtie2 \\
        -x host_index \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --${sensitivity} \\
        -p ${threads} \\
        --un-conc-gz ${sample_id}_unmapped_R%.fastq.gz \\
        -S /dev/null \\
        2> ${sample_id}.bowtie2.log
    """
}

process MULTIQC_HOST_ALIGNMENT {
    tag "MultiQC host alignment"
    label 'process_low'

    conda 'bioconda::multiqc python=3.12'

    input:
    path bowtie2_logs

    output:
    path "host_alignment_multiqc.html", emit: report
    path "host_alignment_multiqc_data", emit: data_dir

    script:
    """
    multiqc . -o . -n host_alignment_multiqc.html
    """
}

// Build the host bowtie2 index once. Call this from the top-level workflow
// (e.g. main_merge.nf) and thread the resulting index channel into every batch
// via HOST_ALIGNMENT, mirroring how silva_classifier_v3v4_ch is shared.
workflow WF_PREPARE_HOST_INDEX {
    take:
    host_fasta  // String: local path or http(s)/ftp URL

    main:
    if (host_fasta == null) {
        error "Host alignment is enabled but params.host_alignment.host_fasta is not set"
    }

    def isUrl = (host_fasta as String) ==~ /^(https?|ftp):\/\/.*/
    if (isUrl) {
        fasta_ch = DOWNLOAD_HOST_FASTA(host_fasta).fasta
    } else {
        fasta_ch = Channel.fromPath(host_fasta, checkIfExists: true)
    }

    // .first() turns the queue emit into a value channel so the index can be
    // broadcast to every batch's HOST_ALIGNMENT call (channels are queues by
    // default and would only deliver to a single consumer).
    index_ch = BUILD_BOWTIE2_INDEX(fasta_ch).index_files.first()

    emit:
    index_files = index_ch
}

// Per-batch alignment + MultiQC. Consumes a prebuilt index channel.
workflow WF_BOWTIE2_HOST_ALIGNMENT {
    take:
    read_pairs_ch    // [sample_id, [r1, r2]]
    host_index_ch    // path channel: emit from WF_PREPARE_HOST_INDEX

    main:
    align_results = ALIGN_HOST_BOWTIE2(read_pairs_ch, host_index_ch)

    multiqc_input_ch = align_results.log.collect()
    multiqc_results = MULTIQC_HOST_ALIGNMENT(multiqc_input_ch)

    // Reshape unmapped pairs into the canonical [sample_id, [r1, r2]] tuple
    // used everywhere else in the pipeline.
    unmapped_reads_ch = align_results.unmapped.map { sample_id, r1, r2 ->
        tuple(sample_id, [r1, r2])
    }

    emit:
    logs           = align_results.log
    unmapped_reads = unmapped_reads_ch
    multiqc_report = multiqc_results.report
    multiqc_data   = multiqc_results.data_dir
}

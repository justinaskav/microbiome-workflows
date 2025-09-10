// Generic cutadapt workflow for trimming sequence reads
// Works with any paired-end sequencing data

process IMPORT_READS_TO_QZA {
    tag "Import ${sampleId}"

    conda params.qiime2_conda_env

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_demux.qza"), emit: demux_artifact

    script:
    """
    # Create directory for CasavaOneEightSingleLanePerSampleDirFmt
    mkdir -p ${sampleId}_casava_import
    
    # Rename reads to match Casava format expected by QIIME2
    cp -L ${reads[0]} ${sampleId}_casava_import/${sampleId}_01_L001_R1_001.fastq.gz
    cp -L ${reads[1]} ${sampleId}_casava_import/${sampleId}_01_L001_R2_001.fastq.gz
    
    # Import as QIIME2 artifact
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path ${sampleId}_casava_import \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path ${sampleId}_demux.qza
    """
}

process TRIM_CUTADAPT {
    tag "Trim ${sampleId}"

    conda params.qiime2_conda_env

    input:
    tuple val(sampleId), path(demuxArtifact), val(barcode_f), val(barcode_r), val(primer_f), val(primer_r)

    output:
    tuple val(sampleId), path("${sampleId}_trimmed.qza"), emit: trimmed_artifact

    script:
    def discard_untrimmed = true
    def minimum_length = 100

    def forw_read_5 = "${barcode_f}${primer_f}"
    def rev_read_5 = "${barcode_r}${primer_r}"

    """
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences ${demuxArtifact} \
      --p-front-f ${forw_read_5} \
      --p-front-r ${rev_read_5} \
      --p-minimum-length ${minimum_length} \
      --p-discard-untrimmed \
      --o-trimmed-sequences ${sampleId}_trimmed.qza \
      --verbose
    """
}

process TRIM_POLYG_SEQUENCES {
    tag "Trim polyG ${sampleId}"

    conda params.qiime2_conda_env

    input:
    tuple val(sampleId), path(trimmedArtifact)

    output:
    tuple val(sampleId), path("${sampleId}_polyg_trimmed.qza"), emit: polyg_trimmed_artifact

    script:
    def minimum_length = 100

    """
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences ${trimmedArtifact} \
      --p-adapter-f GGGGGGGG \
      --p-adapter-r GGGGGGGG \
      --p-minimum-length ${minimum_length} \
      --o-trimmed-sequences ${sampleId}_polyg_trimmed.qza \
      --verbose
    """
}

process EXPORT_QIIME2_SEQUENCES {
    tag "Export ${sampleId}"

    conda params.qiime2_conda_env

    input:
    tuple val(sampleId), path(qza_artifact)

    output:
    tuple val(sampleId), path("${sampleId}_01_L001_R1_001.fastq.gz"), path("${sampleId}_01_L001_R2_001.fastq.gz"), emit: fastq_files

    script:
    """
    # Export the QZA artifact
    qiime tools export --input-path ${qza_artifact} --output-path exported
    
    # Move files to expected names
    for fq in exported/*.fastq.gz; do
        if [[ \$fq == *"forward"* ]] || [[ \$fq == *"_R1_"* ]]; then
            mv "\$fq" "${sampleId}_01_L001_R1_001.fastq.gz"
        elif [[ \$fq == *"reverse"* ]] || [[ \$fq == *"_R2_"* ]]; then
            mv "\$fq" "${sampleId}_01_L001_R2_001.fastq.gz"
        fi
    done
    
    # Cleanup
    rm -rf exported
    """
}

process EXPORT_QIIME2_ARTIFACT {
    tag "Export QIIME2 artifact"

    conda params.qiime2_conda_env

    input:
    path qza_artifact

    output:
    path "exported/*.fastq.gz"

    script:
    """
    # Export the QZA artifact to fastq format
    qiime tools export --input-path ${qza_artifact} --output-path exported
    """
}

process IMPORT_FASTQS_TO_QIIME2 {
    tag "Import FASTQ files to QIIME2"

    conda params.qiime2_conda_env

    input:
    path "*"

    output:
    path "sequences.qza", emit: sequences
    path "sequences_summary.qzv", emit: sequences_summary

    script:
    """
    # Import all files into a single artifact
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path . \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path sequences.qza
    
    # Summarize the sequences
    qiime demux summarize \
      --i-data sequences.qza \
      --o-visualization sequences_summary.qzv
    """
}

process STAGE_FASTQS_FOR_IMPORT {
    tag "Stage fastq files for bulk import"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple path("${sample_id}_01_L001_R1_001.fastq.gz"), path("${sample_id}_01_L001_R2_001.fastq.gz")

    script:
    """
    # Create symbolic links with proper Casava naming
    ln -sf ${reads[0]} ${sample_id}_01_L001_R1_001.fastq.gz
    ln -sf ${reads[1]} ${sample_id}_01_L001_R2_001.fastq.gz
    """
}

process TRIM_MERGED_CUTADAPT {
    tag "Trim merged sequences"

    label "process_medium"
    conda params.qiime2_conda_env

    input:
    path merged_sequences
    val primer_f
    val primer_r

    output:
    path "trimmed_sequences.qza", emit: trimmed_sequences
    path "trimmed_sequences_summary.qzv", emit: trimmed_sequences_summary

    script:
    def minimum_length = 100

    """
    # Apply primer trimming to all sequences at once
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences ${merged_sequences} \
      --p-front-f ${primer_f} \
      --p-front-r ${primer_r} \
      --p-minimum-length ${minimum_length} \
      --p-discard-untrimmed \
      --o-trimmed-sequences primer_trimmed.qza \
      --verbose
    
    # Apply poly-G trimming
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences primer_trimmed.qza \
      --p-adapter-f GGGGGGGGG \
      --p-adapter-r GGGGGGGGG \
      --p-minimum-length ${minimum_length} \
      --o-trimmed-sequences trimmed_sequences.qza \
      --verbose
    
    # Summarize the final trimmed sequences
    qiime demux summarize \
      --i-data trimmed_sequences.qza \
      --o-visualization trimmed_sequences_summary.qzv
    """
}

// Main workflow for cutadapt trimming
workflow WF_CUTADAPT {
    take:
    read_pairs_ch // Channel with sample ID and read pairs
    barcode_primer_info_ch // Channel with barcode/primer info
    trimming_config // Trimming options (skip_polyg_trim, primer_f, primer_r)
    primers

    main:

    if (trimming_config.use_sampleseq) {
        // TRADITIONAL PATH: Individual sample processing for barcoded data
        // This is used when barcode_primer_info is provided

        // Import sequences for each sample into QIIME2 format
        demux_artifacts_ch = IMPORT_READS_TO_QZA(read_pairs_ch)

        // Create a channel with barcode/primer info for each sample
        barcode_primer_info_ch
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                [row.SampleID, row.Barcode_F, row.Barcode_R, row.Primer_F, row.Primer_R]
            }
            .set { sample_info_ch }

        // Prepare demux artifacts for joining with sample info
        // Create channel with both full ID (for downstream) and extracted name (for joining)
        // Extract actual sample name from complex IDs: FDMP20H000840-1a_L1_SimB4 -> SimB4
        demux_for_join_ch = demux_artifacts_ch.map { fullSampleId, artifact ->
            // Extract sample name - handle prefixed patterns but preserve simple names
            def extractedName = fullSampleId.contains('_L1_')
                ? fullSampleId.replaceAll(/.*_L1_/, '')
                : fullSampleId
            [extractedName, fullSampleId, artifact]
        }

        // Smart join with error handling
        // Join on extracted sample name, but preserve full sample ID for downstream processes
        joined_ch = demux_for_join_ch.join(sample_info_ch, by: 0, failOnMismatch: false)

        // Check if join produced results and provide helpful error if not
        joined_count_ch = joined_ch.count()
        joined_count_ch.subscribe { count ->
            if (count == 0) {
                println("ERROR: No sample ID matches found between FASTQ files and SampleSeq data")
                println("This usually means:")
                println("1. Sample names in FASTQ files don't match those in SampleSeq metadata")
                println("2. SampleSeq file is missing or malformed")
                println("3. FASTQ files have unexpected naming pattern")
                throw new Exception("Sample ID join failed - check sample naming consistency")
            }
            else {
                println("Successfully joined ${count} samples between FASTQ and SampleSeq data")
            }
        }

        // Create trimming input channel with proper sample ID handling
        trimming_input_ch = joined_ch.map { extractedName, fullSampleId, artifact, barcode_f, barcode_r, _primer_f, _primer_r ->
            // Use full sample ID for downstream processes (preserves original behavior)
            [
                sampleId: fullSampleId,
                demuxArtifact: artifact,
                barcode_f: barcode_f,
                barcode_r: barcode_r,
                primer_f: _primer_f,
                primer_r: _primer_r,
            ]
        }

        // Run cutadapt for each sample
        trimmed_artifacts_ch = TRIM_CUTADAPT(
            trimming_input_ch
        )

        if (trimming_config.skip_polyg_trim) {
            // If polyG trimming is skipped, just use the trimmed artifacts as they are
            polyg_trimmed_artifacts_ch = trimmed_artifacts_ch
        }
        else {
            // Run polyG trimming on the already trimmed sequences
            polyg_trimmed_artifacts_ch = TRIM_POLYG_SEQUENCES(trimmed_artifacts_ch)
        }

        // Export polyA-trimmed sequences to fastq
        exported_fastqs_ch = EXPORT_QIIME2_SEQUENCES(polyg_trimmed_artifacts_ch)

        // Merge all polyG-trimmed sequences into a single artifact
        merged_results_ch = IMPORT_FASTQS_TO_QIIME2(
            exported_fastqs_ch.map { sample_id, r1, r2 -> [r1, r2] }.flatten().collect()
        )

        trimmed_sequences = merged_results_ch.sequences
        trimmed_sequences_summary = merged_results_ch.sequences_summary
        individual_trimmed = polyg_trimmed_artifacts_ch[1]
    }
    else if (primers.forward != null && primers.reverse != null) {
        // OPTIMIZED PATH: Direct import → cutadapt → DADA2 for same primers
        // This is used when primers are explicitly provided
        staged_fastqs_ch = STAGE_FASTQS_FOR_IMPORT(read_pairs_ch)
            .map { r1, r2 ->
                [r1, r2]
            }
            .flatten()
            .collect()

        // Import all sequences at once
        import_results = IMPORT_FASTQS_TO_QIIME2(staged_fastqs_ch)

        // Apply cutadapt with primers to merged sequences
        trim_results = TRIM_MERGED_CUTADAPT(
            import_results.sequences,
            primers.forward,
            primers.reverse,
        )

        exported_fastqs_ch = EXPORT_QIIME2_ARTIFACT(
            trim_results.trimmed_sequences
        ).flatMap { list_of_files ->
            list_of_files.collect { file ->
                def file_name = file.getName()
                def sample_id = file_name.split('_')[0]
                // Assuming file names are like sample_01_L001_R1_001.fastq.gz
                [sample_id, file]
            }
        }.groupTuple().map { sample_id, files ->
            [sample_id, files[0], files[1]]
        }

        trimmed_sequences = trim_results.trimmed_sequences
        trimmed_sequences_summary = trim_results.trimmed_sequences_summary
        individual_trimmed = Channel.empty()
    }

    emit:
    trimmed_sequences = trimmed_sequences // .qza artifact
    trimmed_sequences_summary = trimmed_sequences_summary //.qzv summary
    individual_trimmed = individual_trimmed // Individual trimmed artifacts
    exported_fastqs = exported_fastqs_ch // Exported fastq files
}

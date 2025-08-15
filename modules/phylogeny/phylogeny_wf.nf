// Workflow for Phylogenetic Tree Construction using QIIME2 align-to-tree-mafft-fasttree pipeline

// Default parameters specific to this workflow module
params.phylogeny_n_threads = 4            // Threads for MAFFT part of the pipeline
params.phylogeny_mask_max_gap_freq = 1.0  // Max gap frequency for masking (1.0 means no filtering by this criterion initially)
params.phylogeny_mask_min_conservation = 0.4// Min conservation for masking (0.4 is a common default)

process ALIGN_AND_TREE_QIIME2 {
    tag "phylogeny_pipeline_${rep_seqs_qza.getBaseName()}"
    conda "${params.qiime2_conda_env}"
    cpus params.phylogeny_n_threads ?: 4

    input:
    path rep_seqs_qza // Representative sequences (QZA)

    output:
    path("alignment.qza"), emit: alignment_qza
    path("masked_alignment.qza"), emit: masked_alignment_qza
    path("unrooted_tree.qza"), emit: unrooted_tree_qza
    path("rooted_tree.qza"), emit: rooted_tree_qza

    script:
    def threads_to_use = task.cpus > 0 ? task.cpus : 1
    // Note: FastTree's own threading is often via OMP_NUM_THREADS,
    // but --p-n-threads in this pipeline is passed to MAFFT.
    // We set OMP_NUM_THREADS as well for FastTree just in case, though the pipeline might manage it.
    """
    export OMP_NUM_THREADS=${threads_to_use}
    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences "${rep_seqs_qza}" \\
        --p-n-threads ${threads_to_use} \\
        --p-mask-max-gap-frequency ${params.phylogeny_mask_max_gap_freq} \\
        --p-mask-min-conservation ${params.phylogeny_mask_min_conservation} \\
        --o-alignment alignment.qza \\
        --o-masked-alignment masked_alignment.qza \\
        --o-tree unrooted_tree.qza \\
        --o-rooted-tree rooted_tree.qza \\
        --verbose
    """
}

workflow WF_PHYLOGENY {
    take:
        rep_seqs_ch // Channel of representative sequences (QZA)

    main:
        ALIGN_AND_TREE_QIIME2 ( rep_seqs_ch )

    emit:
        aligned_sequences = ALIGN_AND_TREE_QIIME2.out.alignment_qza
        masked_alignment = ALIGN_AND_TREE_QIIME2.out.masked_alignment_qza
        unrooted_tree = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
        rooted_tree = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
}

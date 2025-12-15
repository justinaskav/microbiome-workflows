// Phylogeny Subworkflow
// Performs sequence alignment and phylogenetic tree construction

nextflow.enable.dsl = 2

include { ALIGN_AND_TREE_QIIME2 } from '../modules/phylogeny/phylogeny_wf.nf'

workflow PHYLOGENY {
    take:
    rep_seqs  // Representative sequences (.qza)

    main:
    ALIGN_AND_TREE_QIIME2(rep_seqs)

    emit:
    alignment = ALIGN_AND_TREE_QIIME2.out.alignment_qza
    masked_alignment = ALIGN_AND_TREE_QIIME2.out.masked_alignment_qza
    unrooted_tree = ALIGN_AND_TREE_QIIME2.out.unrooted_tree_qza
    rooted_tree = ALIGN_AND_TREE_QIIME2.out.rooted_tree_qza
}

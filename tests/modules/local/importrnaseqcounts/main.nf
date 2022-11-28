#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IMPORTRNASEQCOUNTS } from '../../../../../modules/local/IMPORTRNASEQCOUNTS/main.nf'

workflow test_deseq2_importrnaseqcounts {

    ch_samplesheet = Channel.fromPath("/home-link/iivow01/git/nf-core-modules/modules/nf-core/deseq2/importfeaturecounts/testdata/Sample_preparations.tsv")
    ch_counts = Channel.fromPath("/home-link/iivow01/git/nf-core-modules/modules/nf-core/deseq2/importfeaturecounts/testdata/rnaseq_gene_counts.txt")

    IMPORTRNASEQCOUNTS (
        ch_samplesheet,
        ch_counts
    )
}

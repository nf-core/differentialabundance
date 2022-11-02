#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IMPORTMERGEDCOUNTS } from '../../../../../modules/local/IMPORTMERGEDCOUNTS/main.nf'

workflow test_deseq2_importmergedcounts {

    ch_samplesheet = Channel.fromPath("/home-link/iivow01/git/nf-core-modules/modules/nf-core/deseq2/importfeaturecounts/testdata/Sample_preparations.tsv")
    ch_counts = Channel.fromPath("/home-link/iivow01/git/nf-core-modules/modules/nf-core/deseq2/importfeaturecounts/testdata/merged_gene_counts.txt")

    IMPORTMERGEDCOUNTS (
        ch_samplesheet,
        ch_counts
    )
}

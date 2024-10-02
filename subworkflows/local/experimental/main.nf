//
// Run experimental analysis
//
include { CORRELATION }         from '../correlation/main.nf'
include { DIFFERENTIAL }        from '../differential/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'

workflow EXPERIMENTAL {
    take:
    ch_contrasts
    ch_samplesheet
    ch_counts
    ch_tools

    main:

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    DIFFERENTIAL(
        ch_contrasts,
        ch_samplesheet,
        ch_counts,
        ch_tools
    )
    ch_diff_results = DIFFERENTIAL.out.results
    ch_diff_adjacency = DIFFERENTIAL.out.adjacency

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    CORRELATION(
        ch_counts,
        ch_tools
    )
    ch_corr_matrix = CORRELATION.out.matrix
    ch_corr_adjacency = CORRELATION.out.adjacency

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    ENRICHMENT(
        ch_diff_adjacency,
        ch_corr_adjacency,
        ch_counts
    )
    ch_enriched_cor = ENRICHMENT.out.enriched_cor
    ch_enriched_diff = ENRICHMENT.out.enriched_diff

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: add visualization stuff here

    // do we need to emit anything?
    emit:
    diff_res        = ch_diff_results
    diff_adj        = ch_diff_adjacency
    corr_matrix     = ch_corr_matrix
    corr_adj        = ch_corr_adjacency
    enriched_cor    = ch_enriched_cor
    enriched_cor    = ch_enriched_diff
}

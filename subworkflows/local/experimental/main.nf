//
// Run experimental analysis
//
include { CORRELATION }         from '../correlation/main.nf'
include { DIFFERENTIAL }        from '../differential/main.nf'
include { VARIABLE_SELECTION }  from '../variable_selection/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'


workflow EXPERIMENTAL {
    take:
    ch_contrasts
    ch_samplesheet
    ch_counts
    ch_tools


    main:
    // Perform differential analysis
    DIFFERENTIAL(
        ch_contrasts,
        ch_samplesheet,
        ch_counts,
        ch_tools
    )
    ch_diff_results = DIFFERENTIAL.out.results
    ch_diff_adjacency = DIFFERENTIAL.out.adjacency

    // Perform variable selection
    ch_counts_filtered = VARIABLE_SELECTION(ch_diff_adjacency, ch_counts)

    // Perform correlation analysis
    CORRELATION(
        ch_counts,
        ch_tools,
        ch_counts_filtered
    )
    ch_matrix = CORRELATION.out.matrix
    ch_cor_adjacency = CORRELATION.out.adjacency

    // Perform enrichment analysis
    ENRICHMENT(
        ch_diff_adjacency,
        ch_cor_adjacency,
        ch_counts
    )
    ch_enriched_cor = ENRICHMENT.out.enriched_cor
    ch_enriched_diff = ENRICHMENT.out.enriched_diff

    emit:
    diff_res    = ch_diff_results
    diff_adj    = ch_diff_adjacency
    var_count   = ch_counts_filtered
    corr_matrix = ch_matrix
    corr_adj    = ch_cor_adjacency
    enriched_cor    = ch_enriched_cor
    enriched_cor    = ch_enriched_diff
}

//
// Run experimental analysis
//
include { DIFFERENTIAL }        from '../differential/main.nf'
include { CORRELATION }         from '../correlation/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'

workflow EXPERIMENTAL {
    take:
    ch_contrasts
    ch_samplesheet
    ch_counts
    ch_tools

    main:

    ch_tools.view()

    // initialize empty results channels
    ch_results = Channel.empty()      // differential results - it should be a table
    ch_adjacency = Channel.empty()    // adjacency matrix showing the connections between the genes, with values 1|0
    ch_matrix = Channel.empty()       // correlation matrix
    ch_enriched = Channel.empty()     // output table from enrichment analysis

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    DIFFERENTIAL(
        ch_tools,
        ch_counts,
        ch_samplesheet,
        ch_contrasts
    )
    ch_results = ch_results.mix(DIFFERENTIAL.out.results)
    ch_adjacency = ch_adjacency.mix(DIFFERENTIAL.out.adjacency)

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    CORRELATION(
        ch_tools,
        ch_counts
    )
    ch_matrix = ch_matrix.mix(CORRELATION.out.matrix)
    ch_adjacency = ch_adjacency.mix(CORRELATION.out.adjacency)

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    ENRICHMENT(
        ch_counts,
        ch_results,
        ch_adjacency
    )
    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: call visualization stuff here
}

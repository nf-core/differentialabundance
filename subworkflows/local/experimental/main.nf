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

    // check tools
    ch_tools.view()

    // initialize empty results channels
    ch_results_pairwise = Channel.empty()               // differential results for pairwise analysis - it should be a table
    ch_results_pairwise_filtered = Channel.empty()      // differential results for pairwise analysis - filtered - it should be a table
    ch_results_genewise = Channel.empty()               // differential results for genewise analysis - it should be a table
    ch_results_genewise_filtered = Channel.empty()      // differential results for genewise analysis - filtered - it should be a table
    ch_adjacency = Channel.empty()                      // adjacency matrix showing the connections between the genes, with values 1|0
    ch_matrix = Channel.empty()                         // correlation matrix
    ch_enriched = Channel.empty()                       // output table from enrichment analysis

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    DIFFERENTIAL(
        ch_tools,
        ch_counts,
        ch_samplesheet,
        ch_contrasts
    )
    ch_results_pairwise = ch_results_pairwise.mix(DIFFERENTIAL.out.results_pairwise)
    ch_results_pairwise_filtered = ch_results_pairwise_filtered.mix(DIFFERENTIAL.out.results_pairwise_filtered)
    ch_results_genewise = ch_results_genewise.mix(DIFFERENTIAL.out.results_genewise)
    ch_results_genewise_filtered = ch_results_genewise_filtered.mix(DIFFERENTIAL.out.results_genewise_filtered)
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
        ch_results_genewise,
        ch_results_genewise_filtered,
        ch_adjacency,
        params.gene_sets_files
    )
    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: call visualization stuff here
}

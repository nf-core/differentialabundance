//
// Run experimental analysis
//
include { CORRELATION }         from '../correlation/main.nf'
include { DIFFERENTIAL }        from '../differential/main.nf'
include { VARIABLE_SELECTION }  from '../variable_selection/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'


workflow EXPERIMENTAL {
    take:
    ch_samples_and_matrix // [meta, samplesheet, matrix] que viene de differentialabundance
    ch_tools


    main:
    // Split the ch_samples_and_matrix into one channel for the samplesheet and another for the matrix (PROPD takes them separately).
    ch_samples_and_matrix
    .map {
        meta, samplesheet, counts ->
            [ meta, samplesheet ]
    }
    .set { ch_samplesheet }

    ch_samples_and_matrix
        .map {
            meta, samplesheet, counts ->
                [ meta, counts ]
        }
        .set { ch_counts }

    ch_counts
        .combine(ch_tools)
        .map {
            metacounts, counts, metatools ->
                [ metacounts+metatools, counts ]
        }
        .set { ch_counts_tools }

    // Perform CODA analysis
    ch_out = Channel.empty()

    // Perform differential analysis
    DIFFERENTIAL(ch_counts, ch_tools, ch_samplesheet.collect())
    ch_diff_results = DIFFERENTIAL.out.results
    ch_diff_adjacency = DIFFERENTIAL.out.adjacency

    // Perform variable selection
    ch_counts_filtered = VARIABLE_SELECTION(ch_diff_adjacency, ch_counts)

    // Perform correlation analysis
    CORRELATION(ch_counts, ch_tools, ch_counts_filtered)
    ch_matrix = CORRELATION.out.matrix
    ch_cor_adjacency = CORRELATION.out.adjacency
    ch_out.mix(ch_matrix)

    // Perform enrichment analysis
    ENRICHMENT(ch_diff_adjacency, ch_cor_adjacency, ch_counts)
    ch_enriched_cor = ENRICHMENT.out.enriched_cor
    ch_enriched_diff = ENRICHMENT.out.enriched_diff

    ch_out.mix(ch_enriched_diff, ch_enriched_cor)

    emit:
    output = ch_out
}

//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"

workflow ENRICHMENT {
    take:
    ch_counts
    ch_results_genewise
    ch_results_genewise_filtered
    ch_adjacency
    // TODO: add ch_gm when provided by user, etc.

    main:

    // initialize empty results channels
    ch_enriched = Channel.empty()

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO this should be optional, only run when there is no gene set data provided by user

    MYGENE(ch_counts.take(1))  // only one data is provided to this pipeline
    ch_gmt = MYGENE.out.gmt

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    ch_adjacency
        .filter { it[0]["enr_method"] == "grea" }
        .set { ch_adjacency_grea }

    GREA(ch_adjacency_grea, ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // todo: add gsea here

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // todo: add gprofiler2 here

    emit:
    enriched = ch_enriched
}

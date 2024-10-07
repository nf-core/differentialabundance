//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/nf-core/propr/grea/main.nf"

workflow ENRICHMENT {
    take:
    ch_counts
    ch_results
    ch_adjacency
    // TODO: add ch_gm when provided by user, etc.

    main:

    // initialize empty results channels
    ch_enriched = Channel.empty()

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    // construct the gene set selection
    // TODO this should be optional, only run when there is no gene set data provided by user
    MYGENE(ch_counts)
    ch_gmt = MYGENE.out.gmt

    // GREA method needs adjacency matrix as input
    ch_adjacency
        .filter { it[0]["enr_method"] == "grea" }
        .set { ch_adjacency_grea }

    // run GREA
    GREA(ch_adjacency_grea, ch_gmt.collect())
    ch_enriched = ch_enriched.mix(GREA.out.enrichedGO)

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // todo: add gsea here
    // then we need to add the corresponding input channels to this subworkflow

    emit:
    enriched = ch_enriched
}

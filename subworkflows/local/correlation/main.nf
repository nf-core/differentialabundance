//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/nf-core/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_counts
    ch_tools

    main:

    // initialize empty results channels
    ch_results   = Channel.empty()
    ch_adjacency = Channel.empty()

    // branch tools to select the correct correlation analysis method
    ch_counts
        .combine(ch_tools)
        .map {
            metacounts, counts, metatools ->
                [ metacounts+metatools, counts ]
        }
        .branch {
            propr:  it[0]["cor_method"] == "propr"
        }
        .set { ch_counts_cor }

    // ----------------------------------------------------
    // Perform correlation analysis with propr
    // ----------------------------------------------------

    PROPR(ch_counts_cor.propr)
    ch_matrix = PROPR.out.matrix
    ch_adjacency = PROPR.out.adj

    emit:
    matrix = ch_matrix
    adjacency = ch_adjacency
}

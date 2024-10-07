//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/nf-core/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_tools
    ch_counts

    main:

    // initialize empty results channels
    ch_results   = Channel.empty()
    ch_adjacency = Channel.empty()

    // branch tools to select the correct correlation analysis method
    ch_tools
        .branch {
            propr:  it[0]["cor_method"] == "propr"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform correlation analysis with propr
    // ----------------------------------------------------

    ch_counts
        .combine(ch_tools_single.propr)
        .map {
            metacounts, counts, metatools ->
                [ metacounts+metatools, counts ]
        }
        .set { ch_counts_propr }

    PROPR(ch_counts_propr)
    ch_matrix = PROPR.out.matrix
    ch_adjacency = PROPR.out.adj

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    matrix = ch_matrix
    adjacency = ch_adjacency
}

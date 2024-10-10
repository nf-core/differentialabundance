//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/local/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_tools
    ch_counts

    main:

    // initialize empty results channels
    ch_matrix   = Channel.empty()
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
    ch_matrix = ch_matrix.mix(PROPR.out.matrix)
    ch_adjacency = ch_adjacency.mix(PROPR.out.adjacency)

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    matrix = ch_matrix
    adjacency = ch_adjacency
}

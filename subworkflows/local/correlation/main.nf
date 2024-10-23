//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/local/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_tools        // [ correlation_map ] with the keys: cor_method, args_cor
    ch_counts

    main:

    // initialize empty results channels
    ch_matrix   = Channel.empty()
    ch_adjacency = Channel.empty()

    // branch tools to select the correct correlation analysis method
    ch_tools
        .branch {
            propr:  it["cor_method"] == "propr"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform correlation analysis with propr
    // ----------------------------------------------------

    ch_counts
        .combine(ch_tools_single.propr)
        .map {
            metacounts, counts, metatools ->
                input:   [ metacounts+metatools, counts ]
        }
        .set { ch_counts_propr }

    PROPR(ch_counts_propr.unique())
    ch_matrix    = PROPR.out.matrix.mix(ch_matrix)
    ch_adjacency = PROPR.out.adjacency.mix(ch_adjacency)

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    matrix    = ch_matrix
    adjacency = ch_adjacency
}

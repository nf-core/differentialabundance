//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/local/propr/propr/main.nf"

def correct_meta_data = { meta, data, pathway ->
    def meta_clone = meta.clone() + pathway
    meta_clone.remove('cor_method')
    meta_clone.remove('args_cor')
    return [meta_clone, data]
}

workflow CORRELATION {
    take:
    ch_tools        // [ pathway_name, correlation_map ]
    ch_counts

    main:

    // initialize empty results channels
    ch_matrix   = Channel.empty()
    ch_adjacency = Channel.empty()

    // branch tools to select the correct correlation analysis method
    ch_tools
        .branch {
            propr:  it[1]["cor_method"] == "propr"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform correlation analysis with propr
    // ----------------------------------------------------

    ch_counts
        .combine(ch_tools_single.propr)
        .multiMap {
            metacounts, counts, pathway, metatools ->
                input:   [ metacounts+metatools, counts ]
                pathway: [ metacounts+metatools, pathway ]
        }
        .set { ch_counts_propr }

    PROPR(ch_counts_propr.input.unique())
    ch_matrix    = PROPR.out.matrix
                        .join(ch_counts_propr.pathway).map(correct_meta_data).mix(ch_matrix)
    ch_adjacency = PROPR.out.adjacency
                        .join(ch_counts_propr.pathway).map(correct_meta_data).mix(ch_adjacency)

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    matrix    = ch_matrix
    adjacency = ch_adjacency
}

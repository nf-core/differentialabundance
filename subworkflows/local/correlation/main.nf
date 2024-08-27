//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/nf-core/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_counts
    ch_tools
    ch_counts_filtered

    main:
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

    // Create a branch of the channel to retrieve the normal counts when there is no variable selection.
    ch_counts_cor.propr
        .branch{
            no_sel: it[0]["sel_method"] == null
            sel:    it[0]["sel_method"] != null
        }
        .set { ch_counts_selection }

    ch_propr = ch_counts_filtered.mix(ch_counts_selection.no_sel)

    PROPR(ch_propr)
    ch_matrix = PROPR.out.matrix
    ch_adjacency = PROPR.out.adj

    emit:
    matrix = ch_matrix
    adjacency = ch_adjacency

}

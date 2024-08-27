//
// Perform variable selection
//
include { FILTERVAR } from "../../../modules/local/filtervar/main.nf"

workflow VARIABLE_SELECTION {
    take:
    ch_adj      //meta_tools, adj
    ch_counts   //meta_id, counts

    main:
    ch_counts
        .map {
            metacounts, counts ->
                [counts]
        }
        .combine(ch_adj)
        .map{
            counts, meta, adj ->
                [ meta, counts, adj]
        }
        .branch {
            filtervar:    it[0]["sel_method"] == "filtervar"
            deseqfilter:  it[0]["sel_method"] == "deseqfilter"
        }
        .set { ch_counts_adj_sel }

    FILTERVAR(ch_counts_adj_sel.filtervar)

    ch_counts_cor = FILTERVAR.out.count

    emit:
    count = ch_counts_cor
}

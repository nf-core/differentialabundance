// include modules
include { FILTERVAR } from "../../modules/local/filtervar/main.nf"

workflow VARIABLE_SELECTION {
    take:
    ch_adj//meta_tools, adj
    ch_counts //meta_id, counts


    main:
    ch_counts
        .map {
            metacounts, counts ->
                [counts]
        }
        .combine(ch_adj)
        //.view()
        .map{
            counts, meta, adj -> 
                [ meta, counts, adj]
        }
        //.view()
        .branch {
            filtervar:    it[0]["sel_method"] == "filtervar"
            deseqfilter:  it[0]["sel_method"] == "deseqfilter"
        }
        .set { ch_counts_adj_sel }

    //ch_counts_adj_sel.nofilter.view()


    FILTERVAR(ch_counts_adj_sel.filtervar)

    ch_counts_cor = FILTERVAR.out.count
    //ch_counts_cor.view()

    
    emit:
    count = ch_counts_cor
}
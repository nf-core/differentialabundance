// include modules
include {PROPR_PROPD as PROPD} from "../../modules/nf-core/propr/propd/main.nf"


workflow DIFFERENTIAL {
    take:
    ch_counts
    ch_tools
    ch_samplesheet

    main:
    ch_counts
        .combine(ch_tools)
        .map {
            metacounts, counts, meta ->
                [ metacounts+meta, counts ]
        }
        //.view()
        .branch {
            propd:    it[0]["diff_method"] == "propd"
            deseq2: it[0]["diff_method"] == "deseq2"
        }
        .set { ch_counts_tools }
     
    PROPD(ch_counts_tools.propd, ch_samplesheet)
    ch_results = PROPD.out.results
    ch_adjacency = PROPD.out.adj

    emit:
    results = ch_results
    adjacency = ch_adjacency

}
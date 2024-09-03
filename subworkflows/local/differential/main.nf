//
// Perform differential analysis
//
include {PROPR_PROPD as PROPD} from "../../../modules/nf-core/propr/propd/main.nf"

workflow DIFFERENTIAL {
    take:
    ch_contrasts    // [meta, contrast_variable, reference, target]
    ch_counts
    ch_tools
    ch_samplesheet

    main:
    ch_counts
        .combine(ch_tools)
        .combine(ch_contrasts)
        .map {
            meta_counts, counts, tools, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_counts.clone() + tools.clone()
                meta.args_diff = meta.args_diff + " --group_col $contrast_variable"
                [ meta, counts ]
        }
        .branch {
            propd:  it[0]["diff_method"] == "propd"
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

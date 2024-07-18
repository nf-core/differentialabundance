// include modules

include { PROPR_GREA as GREA_DIFF } from "../../modules/nf-core/propr/grea/main.nf"
include { PROPR_GREA as GREA_COR } from "../../modules/nf-core/propr/grea/main.nf"
include { MYGENE } from "../../modules/nf-core/mygene/main.nf"


workflow ENRICHMENT {
    take:
    ch_diff_adjacency
    ch_cor_adjacency 
    ch_counts 


    main:

    MYGENE(ch_counts)
    ch_gmt = MYGENE.out.gmt


    ch_diff_adjacency
        .branch {
            grea: it[0]["enr_diff_method"] == "grea"
            gsea: it[0]["enr_diff_method"] == "gsea"
        }
        .set { ch_diff_grea }

    GREA_DIFF(ch_diff_grea.grea, ch_gmt.collect())
    ch_enriched_diff = GREA_DIFF.out.enrichedGO

    ch_cor_adjacency
        .branch {
            grea: it[0]["enr_cor_method"] == "grea"        }
        .set { ch_cor_grea }

    ch_cor_grea.grea.view()
    ch_diff_grea.grea.view()

    GREA_COR(ch_cor_grea.grea, ch_gmt.collect())
    ch_enriched_cor = GREA_COR.out.enrichedGO


    emit:
    enriched_diff = ch_enriched_diff
    enriched_cor = ch_enriched_cor
}
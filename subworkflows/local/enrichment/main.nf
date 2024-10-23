//
// Perform enrichment analysis
//
include { MYGENE } from "../../../modules/nf-core/mygene/main.nf"
include { PROPR_GREA as GREA } from "../../../modules/local/propr/grea/main.nf"
include { GPROFILER2_GOST } from "../../../modules/nf-core/gprofiler2/gost/main.nf"

workflow ENRICHMENT {
    take:
    ch_tools        // [ pathway_name, enrichment_map ]
    ch_counts
    ch_results_genewise
    ch_results_genewise_filtered
    ch_adjacency
    gene_sets_file

    main:

    // initialize empty results channels
    ch_enriched = Channel.empty()

    // ----------------------------------------------------
    // Construct gene set database
    // ----------------------------------------------------

    // TODO for the moment, it only works with one gene sets file
    // we may want to support multiple gene sets files

    // TODO either we need to make the mygene module more complete to parse GMT file from different databases,
    // or use other alternatives like for example the gene set parsed from Gprofiler2. (The problem with this is that it only works for certain organisms)

    if (gene_sets_file == null) {
        MYGENE(ch_counts.take(1))
        ch_gene_sets = MYGENE.out.gmt
    } else {
        Channel.of(gene_sets_file)
            .map { file(it, checkIfExists: true) }
            .map { it -> [[id: it.name], it]}
            .set { ch_gene_sets }
    }

    // TODO the gene ids in the gene sets should be the same as the gene ids in the counts table
    // otherwise GREA won't work properly
    // We may want to add a step to check this and convert between gene_id/gene_name if necessary

    // ----------------------------------------------------
    // Perform enrichment analysis with GREA
    // ----------------------------------------------------

    ch_adjacency
        .map { meta, matrix -> [meta.subMap(["pathway_name"]), meta, matrix] }
        .join(ch_tools, by: [0])
        .map {
            pathway_name, meta, matrix, meta_tools ->
                def new_meta = meta.clone() + meta_tools.clone()
                [ new_meta, matrix ]
            }
        .branch {
            grea:  it[0]["enr_method"] == "grea"
            gsea: it[0]["enr_method"] == "gsea"
        }
        .set { ch_adjacency }

    GREA(ch_adjacency.grea, ch_gene_sets.collect())
    ch_enriched = ch_enriched.mix(GREA.out.results)

    // ----------------------------------------------------
    // Perform enrichment analysis with gprofiler2
    // ----------------------------------------------------

    // TODO here we assume that ch_gene_sets and ch_counts have only one element, so they can be combined interchangeably
    // but in the future, with differently processed data, etc., we may want to combine them based on some criteria
    ch_results_genewise_filtered
        .filter { it[0]["enr_method"] == "gprofiler2" }
        .combine(ch_gene_sets)
        .combine(ch_counts)
        .multiMap { meta_results, results, meta_gene_sets, gene_sets, meta_counts, counts ->
            de : [meta_results, results]
            gmt : [gene_sets]
            background : [counts]
        }
        .set{ ch_enrichment_gprofiler2 }

    GPROFILER2_GOST(
        ch_enrichment_gprofiler2.de,
        ch_enrichment_gprofiler2.gmt,
        ch_enrichment_gprofiler2.background
    )

    // ----------------------------------------------------
    // Perform enrichment analysis with GSEA
    // ----------------------------------------------------

    // todo: add gsea here

    emit:
    enriched = ch_enriched
}

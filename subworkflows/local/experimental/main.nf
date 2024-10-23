//
// Run experimental analysis
//
include { DIFFERENTIAL }        from '../differential/main.nf'
include { CORRELATION }         from '../correlation/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'

def postprocess_subworkflow_output( ch_results, ch_tools_args ) {
    // 1) join results with pathway data
    // 2) clean up meta data by removing tool arguments and adding pathway name
    return ch_results
        .combine(ch_tools_args)
        .filter{ meta, data, pathway, arg_map -> meta.subMap(arg_map.keySet()) == arg_map }
        .map{ meta, data, pathway, arg_map ->
            def meta_clone = meta.clone() + pathway
            meta_clone.removeAll{it.key in arg_map.keySet()}
            return [meta_clone, data]
        }
}

workflow EXPERIMENTAL {
    take:
    ch_contrasts    // [ meta, contrast_variable, reference, target ]
    ch_samplesheet  // [ meta, samplesheet ]
    ch_counts       // [ meta, counts]
    ch_tools        // [ pathway_name, differential_map, correlation_map, enrichment_map ]

    main:

    // split toolsheet into channels
    ch_tools
        .multiMap{
            pathway_name, differential_map, correlation_map, enrichment_map ->
                diff: [ pathway_name, differential_map ]
                corr: [ pathway_name, correlation_map ]
                enr:  [ pathway_name, enrichment_map ]
        }
        .set{ ch_tools }


    // initialize empty results channels
    ch_results_pairwise = Channel.empty()               // differential results for pairwise analysis - it should be a table
    ch_results_pairwise_filtered = Channel.empty()      // differential results for pairwise analysis - filtered - it should be a table
    ch_results_genewise = Channel.empty()               // differential results for genewise analysis - it should be a table
    ch_results_genewise_filtered = Channel.empty()      // differential results for genewise analysis - filtered - it should be a table
    ch_adjacency = Channel.empty()                      // adjacency matrix showing the connections between the genes, with values 1|0
    ch_matrix = Channel.empty()                         // correlation matrix
    ch_enriched = Channel.empty()                       // output table from enrichment analysis

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    DIFFERENTIAL(
        ch_tools.diff.map{ it[1] },
        ch_counts,
        ch_samplesheet,
        ch_contrasts
    )
    ch_results_pairwise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise,ch_tools.diff).mix(ch_results_pairwise)
    ch_results_pairwise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise_filtered,ch_tools.diff).mix(ch_results_pairwise_filtered)
    ch_results_genewise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise,ch_tools.diff).mix(ch_results_genewise)
    ch_results_genewise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise_filtered,ch_tools.diff).mix(ch_results_genewise_filtered)
    ch_adjacency                 = postprocess_subworkflow_output(DIFFERENTIAL.out.adjacency,ch_tools.diff).mix(ch_adjacency)

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    CORRELATION(
        ch_tools.corr.map{ it[1] },
        ch_counts
    )
    ch_matrix    = postprocess_subworkflow_output(CORRELATION.out.matrix,ch_tools.corr).mix(ch_matrix)
    ch_adjacency = postprocess_subworkflow_output(CORRELATION.out.adjacency,ch_tools.corr).mix(ch_adjacency)

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    ENRICHMENT(
        ch_tools.enr,
        ch_counts,
        ch_results_genewise,
        ch_results_genewise_filtered,
        ch_adjacency
    )
    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: call visualization stuff here
}

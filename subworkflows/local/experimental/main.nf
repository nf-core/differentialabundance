//
// Run experimental analysis
//
include { DIFFERENTIAL }        from '../differential/main.nf'
include { CORRELATION }         from '../correlation/main.nf'
include { ENRICHMENT }          from '../enrichment/main.nf'

def preprocess_subworkflow_output( ch_input, ch_tools_args, method_field_name) {
    // add method arguments to channel meta
    return ch_input
        .combine(ch_tools_args)
        .filter{ meta, input, pathway, arg_map -> meta["pathway_name"] ? meta["pathway_name"] == pathway["pathway_name"] : true }
        .map{ meta, input, pathway, arg_map ->
            def meta_clone = meta.clone() + pathway + arg_map.clone()
            def method = meta_clone.remove(method_field_name)
            return [meta_clone, input, method]
        }
}

def postprocess_subworkflow_output( ch_results, field_name ) {
    // clean up meta data by removing tool arguments
    return ch_results
        .map{ meta, data ->
            def meta_clone = meta.clone()
            meta_clone.removeAll{it.key in field_name}
            return [meta_clone, data]
        }
}

workflow EXPERIMENTAL {
    take:
    ch_contrasts    // [ meta, contrast_variable, reference, target ]
    ch_samplesheet  // [ meta, samplesheet ]
    ch_featuresheet // [ meta, featuresheet ]
    ch_gene_sets
    ch_counts       // [ meta, counts]
    ch_tools        // [ pathway_name, differential_map, correlation_map, enrichment_map ]
    ch_transcript_lengths
    ch_control_features

    main:

    ch_versions = Channel.empty()

    // split toolsheet into channels
    ch_tools.count()
        .combine(ch_tools)
        .multiMap{
            n, pathway, differential_map, correlation_map, enrichment_map ->
                def pathway_name = n == 1 ? ["pathway_name":""] : pathway
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
    ch_matrix    = Channel.empty()                      // correlation matrix
    ch_enriched  = Channel.empty()                      // output table from enrichment analysis

    // ----------------------------------------------------
    // DIFFERENTIAL ANALYSIS BLOCK
    // ----------------------------------------------------

    preprocess_subworkflow_output(ch_counts, ch_tools.diff, "diff_method")
        .set{ ch_counts_diff }
    DIFFERENTIAL(
        ch_counts_diff,
        ch_samplesheet,
        ch_contrasts,
        ch_transcript_lengths,
        ch_control_features
    )
    ch_results_pairwise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise,["method", "args_diff"]).mix(ch_results_pairwise)
    ch_results_pairwise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_pairwise_filtered,["method", "args_diff"]).mix(ch_results_pairwise_filtered)
    ch_results_genewise          = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise,["method", "args_diff"]).mix(ch_results_genewise)
    ch_results_genewise_filtered = postprocess_subworkflow_output(DIFFERENTIAL.out.results_genewise_filtered,["method", "args_diff"]).mix(ch_results_genewise_filtered)
    ch_adjacency                 = postprocess_subworkflow_output(DIFFERENTIAL.out.adjacency,["method", "args_diff"]).mix(ch_adjacency)
    ch_versions                  = ch_versions.mix(DIFFERENTIAL.out.versions)

    // ----------------------------------------------------
    // CORRELATION ANALYSIS BLOCK
    // ----------------------------------------------------

    preprocess_subworkflow_output(ch_counts, ch_tools.corr, "cor_method")
        .set{ ch_counts_corr }

    CORRELATION(
        ch_counts_corr
    )
    ch_matrix    = postprocess_subworkflow_output(CORRELATION.out.matrix,["method", "args_cor"]).mix(ch_matrix)
    ch_adjacency = postprocess_subworkflow_output(CORRELATION.out.adjacency,["method", "args_cor"]).mix(ch_adjacency)
    ch_versions  = ch_versions.mix(CORRELATION.out.versions)

    // ----------------------------------------------------
    // FUNCTIONAL ENRICHMENT BLOCK
    // ----------------------------------------------------

    preprocess_subworkflow_output(ch_counts, ch_tools.enr, "enr_method")
        .set{ ch_counts_enr }
    preprocess_subworkflow_output(ch_results_genewise, ch_tools.enr, "enr_method")
        .set{ ch_results_genewise_enr }
    preprocess_subworkflow_output(ch_results_genewise_filtered, ch_tools.enr, "enr_method")
        .set{ ch_results_genewise_filtered_enr }
    preprocess_subworkflow_output(ch_adjacency, ch_tools.enr, "enr_method")
        .set{ ch_adjacency_enr }

    ENRICHMENT(
        ch_counts_enr,
        ch_results_genewise_enr,
        ch_results_genewise_filtered_enr,
        ch_adjacency_enr,
        ch_contrasts,
        ch_samplesheet,
        ch_featuresheet,
        ch_gene_sets
    )
    ch_enriched = ch_enriched.mix(ENRICHMENT.out.enriched)
    ch_versions = ch_versions.mix(ENRICHMENT.out.versions)

    // ----------------------------------------------------
    // VISUALIZATION BLOCK
    // ----------------------------------------------------

    // TODO: call visualization stuff here

    emit:
    versions = ch_versions
}

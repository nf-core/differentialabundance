//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL  } from '../../../modules/nf-core/deseq2/differential/main'

def correct_meta_data = { meta, data, pathway ->
    def meta_clone = meta.clone() + pathway
    meta_clone.remove('diff_method')
    meta_clone.remove('args_diff')
    return [meta_clone, data]
}

workflow DIFFERENTIAL {
    take:
    ch_tools        // [ pathway_name, differential_map ]
    ch_counts
    ch_samplesheet
    ch_contrasts    // [meta, contrast_variable, reference, target]

    main:

    // initialize empty results channels
    ch_results_pairwise          = Channel.empty()
    ch_results_pairwise_filtered = Channel.empty()
    ch_results_genewise          = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_adjacency                 = Channel.empty()

    // branch tools to select the correct differential analysis method
    ch_tools
        .branch {
            propd:  it[1]["diff_method"] == "propd"
            deseq2: it[1]["diff_method"] == "deseq2"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform differential analysis with propd
    // ----------------------------------------------------

    ch_counts
        .join(ch_samplesheet)
        .combine(ch_contrasts)
        .combine(ch_tools_single.propd)
        .multiMap {
            meta_data, counts, samplesheet, meta_contrast, contrast_variable, reference, target, pathway, meta_tools ->
                def meta = meta_data.clone() + ['contrast': meta_contrast.id] + meta_tools.clone()
                input:   [ meta, counts, samplesheet, contrast_variable, reference, target ]
                pathway: [ meta, pathway ]
        }
        .set { ch_propd }

    PROPD(ch_propd.input.unique())
    ch_results_pairwise          = PROPD.out.results
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_pairwise)
    ch_results_pairwise_filtered = PROPD.out.results_filtered
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_pairwise_filtered)
    ch_results_genewise          = PROPD.out.connectivity
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_genewise)
    ch_results_genewise_filtered = PROPD.out.hub_genes
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_results_genewise_filtered)
    ch_adjacency                 = PROPD.out.adjacency
                                        .join(ch_propd.pathway).map(correct_meta_data).mix(ch_adjacency)

    // ----------------------------------------------------
    // Perform differential analysis with DESeq2
    // ----------------------------------------------------

    // ToDo: In order to use deseq2 the downstream processes need to be updated to process the output correctly
    // if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = [[],[]] }
    // if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = [[],[]] }

    // ch_samplesheet
    //     .join(ch_counts)
    //     .first()
    //     .combine(ch_tools_single.deseq2)
    //     .set { ch_counts_deseq2 }

    // DESEQ2_DIFFERENTIAL (
    //     ch_contrasts,
    //     ch_counts_deseq2,
    //     ch_control_features,
    //     ch_transcript_lengths
    // )
    // ch_results = ch_results
    //     .mix(DESEQ2_DIFFERENTIAL.out.results)

    emit:
    results_pairwise          = ch_results_pairwise
    results_pairwise_filtered = ch_results_pairwise_filtered
    results_genewise          = ch_results_genewise
    results_genewise_filtered = ch_results_genewise_filtered
    adjacency                 = ch_adjacency
}

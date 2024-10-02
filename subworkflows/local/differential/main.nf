//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/nf-core/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL  } from '../../../modules/nf-core/deseq2/differential/main'


workflow DIFFERENTIAL {
    take:
    ch_contrasts    // [meta, contrast_variable, reference, target]
    ch_samplesheet
    ch_counts
    ch_tools

    main:

    // initialize empty results channels
    ch_results   = Channel.empty()
    ch_adjacency = Channel.empty()

    // branch tools to select the correct differential analysis method
    ch_tools
        .branch {
            propd:  it[0]["diff_method"] == "propd"
            deseq2: it[0]["diff_method"] == "deseq2"
        }
        .set { ch_tools_single }

    // ----------------------------------------------------
    // Perform differential analysis with propd
    // ----------------------------------------------------

    ch_counts
        .combine(ch_tools_single.propd)
        .combine(ch_contrasts)
        .map {
            meta_counts, counts, tools, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_counts.clone() + tools.clone()
                meta.args_diff = (meta.args_diff ?: "") + " --group_col $contrast_variable"  // TODO parse the toolsheet with the ext.arg from modules.config at the beginning of the experimental workflow
                [ meta, counts ]
        }
        .unique()
        .set { ch_counts_propd }

    PROPD(
        ch_counts_propd,
        ch_samplesheet.first()
    )
    ch_results = ch_results.mix(PROPD.out.results)
    ch_adjacency = ch_adjacency.mix(PROPD.out.adj)

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
    results   = ch_results
    adjacency = ch_adjacency
}

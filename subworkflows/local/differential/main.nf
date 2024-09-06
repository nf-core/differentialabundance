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
    ch_counts
        .join(ch_samplesheet)
        .first()
        .combine(ch_tools)
        .combine(ch_contrasts)
        .map {
            meta_counts, counts, samplesheet, tools, meta_contrast, contrast_variable, reference, target ->
                def meta = meta_counts.clone() + tools.clone()
                meta.args_diff = (meta.args_diff ?: "") + " --group_col $contrast_variable"
                [ meta, samplesheet, counts ]
        }
        .unique()
        .branch {
            propd:  it[0]["diff_method"] == "propd"
                return [it[0], it[2]]
            deseq2: it[0]["diff_method"] == "deseq2"
        }
        .set { ch_counts_tools }

    PROPD(
        ch_counts_tools.propd,
        ch_samplesheet.first()
    )
    ch_results   = PROPD.out.results
    ch_adjacency = PROPD.out.adj

    // ToDo: In order to use deseq2 the downstream processes need to be updated to process the output correctly
    // if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = [[],[]] }
    // if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = [[],[]] }

    // DESEQ2_DIFFERENTIAL (
    //     ch_contrasts,
    //     ch_counts_tools.deseq2,
    //     ch_control_features,
    //     ch_transcript_lengths
    // )
    // ch_results = ch_results
    //     .mix(DESEQ2_DIFFERENTIAL.out.results)

    emit:
    results   = ch_results
    adjacency = ch_adjacency

}

//
// Perform correlation analysis
//
include {PROPR_PROPR as PROPR} from "../../../modules/local/propr/propr/main.nf"

workflow CORRELATION {
    take:
    ch_counts       // [[meta_input], counts, analysis method]

    main:

    // initialize empty results channels
    ch_matrix    = Channel.empty()
    ch_adjacency = Channel.empty()
    ch_versions  = Channel.empty()

    // branch tools to select the correct correlation analysis method
    inputs = ch_counts
        .map { meta, abundance, analysis_method ->
            [ meta + ['method': analysis_method ], abundance]
        }

    // ----------------------------------------------------
    // Perform correlation analysis with propr
    // ----------------------------------------------------
    PROPR( inputs.filter{it[0].method == 'propr'})
    ch_matrix    = PROPR.out.matrix.mix(ch_matrix)
    ch_adjacency = PROPR.out.adjacency.mix(ch_adjacency)
    ch_versions  = ch_versions.mix(PROPR.out.versions)

    // TODO: divide propr module into cor, propr, pcor, pcorbshrink, etc.

    emit:
    matrix    = ch_matrix     // channel: [ csv ]
    adjacency = ch_adjacency  // channel: [ csv ]
    versions  = ch_versions   // channel: [ versions.yml ]
}

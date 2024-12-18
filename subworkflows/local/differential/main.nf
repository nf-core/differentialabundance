//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential/main'
include { DESEQ2_DIFFERENTIAL as DESEQ2_NORM } from "../../../modules/nf-core/deseq2/differential/main"
include { LIMMA_DIFFERENTIAL } from '../../../modules/nf-core/limma/differential/main'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_LIMMA } from '../../../modules/local/filter_difftable'
include { FILTER_DIFFTABLE as FILTER_DIFFTABLE_DESEQ2 } from '../../../modules/local/filter_difftable'

workflow DIFFERENTIAL {
    take:
    ch_input              // [[meta_input], counts, analysis method]

    ch_samplesheet        // [ meta_exp, samplesheet ]
    ch_contrasts          // [ meta_contrast, contrast_variable, reference, target ]
    ch_transcript_lengths // [ meta_exp, transcript_lengths]
    ch_control_features   // [ meta_exp, control_features]

    main:

    // initialize empty results channels
    ch_results_pairwise          = Channel.empty()
    ch_results_pairwise_filtered = Channel.empty()
    ch_results_genewise          = Channel.empty()
    ch_results_genewise_filtered = Channel.empty()
    ch_adjacency                 = Channel.empty()
    ch_versions                  = Channel.empty()

    // We need to cross the things we're iterating
    inputs = ch_input
        .combine(ch_samplesheet)
        .filter{ meta_input, abundance, analysis_method, meta_exp, samplesheet -> meta_input.subMap(meta_exp.keySet()) == meta_exp }
        .combine(ch_contrasts)
        .multiMap { meta_input, abundance, analysis_method, meta_exp, samplesheet, meta_contrast, variable, reference, target ->
            def meta = meta_input.clone() + meta_contrast.clone() + ['method': analysis_method ]                                                    // ToDo: If modules are updated to combine their input metas, then this merge could be removed
            samples_and_matrix: [ meta, samplesheet, abundance]
            contrasts: [ meta, variable, reference, target]
            propd_input: [ meta + ['contrast': meta_contrast.id, 'method': analysis_method], abundance, samplesheet, variable, reference, target ]  // ToDo: remove this, once propd is updated and uses standardized inputs
        }

    // ----------------------------------------------------
    // Perform differential analysis with propd
    // ----------------------------------------------------

    // TODO propd currently don't support blocking, so we should not run propd with same contrast_variable, reference and target,
    // but different blocking variable, since it will simply run the same analysis again.

    PROPD( inputs.propd_input.filter{it[0].method == 'propd'})
    ch_results_pairwise          = PROPD.out.results.mix(ch_results_pairwise)
    ch_results_pairwise_filtered = PROPD.out.results_filtered.mix(ch_results_pairwise_filtered)
    ch_results_genewise          = PROPD.out.connectivity.mix(ch_results_genewise)
    ch_results_genewise_filtered = PROPD.out.hub_genes.mix(ch_results_genewise_filtered)
    ch_adjacency                 = PROPD.out.adjacency.mix(ch_adjacency)
    ch_versions                  = PROPD.out.versions.mix(ch_versions)

    // ----------------------------------------------------
    // Perform differential analysis with DESeq2
    // ----------------------------------------------------

    // do we need this process DESEQ2_NORM?
    DESEQ2_NORM (
            inputs.contrasts.filter{it[0].method == 'deseq2'}.first(),
            inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
            ch_control_features,
            ch_transcript_lengths
        )

    DESEQ2_DIFFERENTIAL (
            inputs.contrasts.filter{it[0].method == 'deseq2'},
            inputs.samples_and_matrix.filter{it[0].method == 'deseq2'},
            ch_control_features,
            ch_transcript_lengths
        )

    ch_norm_deseq2         = DESEQ2_NORM.out.normalised_counts
    ch_differential_deseq2 = DESEQ2_DIFFERENTIAL.out.results
    ch_model_deseq2        = DESEQ2_DIFFERENTIAL.out.model
    ch_versions            = DESEQ2_DIFFERENTIAL.out.versions.mix(ch_versions)

    // TODO modify the module to accept these parameters as meta/ext.args in the same way how propd does
    ch_logfc_deseq2 = Channel.value([ "log2FoldChange", params.differential_min_fold_change ])
    ch_padj_deseq2 = Channel.value([ "padj", params.differential_max_qval ])

    FILTER_DIFFTABLE_DESEQ2(
        ch_differential_deseq2,
        ch_logfc_deseq2,
        ch_padj_deseq2
    )

    ch_results_genewise          = DESEQ2_DIFFERENTIAL.out.results.mix(ch_results_genewise)
    ch_results_genewise_filtered = FILTER_DIFFTABLE_DESEQ2.out.filtered.mix(ch_results_genewise_filtered)
    ch_versions                  = FILTER_DIFFTABLE_DESEQ2.out.versions.mix(ch_versions)

    // ----------------------------------------------------
    // Perform differential analysis with limma
    // ----------------------------------------------------

    // run limma
    LIMMA_DIFFERENTIAL(
        inputs.contrasts.filter{it[0].method == 'limma'},
        inputs.samples_and_matrix.filter { it[0].method == 'limma' }
    )
    ch_versions = LIMMA_DIFFERENTIAL.out.versions.mix(ch_versions)

    // filter results
    // note that these are column names specific for limma output table
    // TODO modify the module to accept these parameters as meta/ext.args in the same way how propd does
    ch_logfc_limma = Channel.value([ "logFC", params.differential_min_fold_change ])
    ch_padj_limma = Channel.value([ "adj.P.Val", params.differential_max_qval ])
    FILTER_DIFFTABLE_LIMMA(
        LIMMA_DIFFERENTIAL.out.results,
        ch_logfc_limma,
        ch_padj_limma
    )

    // collect results
    ch_results_genewise          = LIMMA_DIFFERENTIAL.out.results.mix(ch_results_genewise)
    ch_results_genewise_filtered = FILTER_DIFFTABLE_LIMMA.out.filtered.mix(ch_results_genewise_filtered)
    ch_versions                  = FILTER_DIFFTABLE_LIMMA.out.versions.mix(ch_versions)

    emit:
    results_pairwise          = ch_results_pairwise           // channel: [ tsv ]
    results_pairwise_filtered = ch_results_pairwise_filtered  // channel: [ tsv ]
    results_genewise          = ch_results_genewise           // channel: [ tsv ]
    results_genewise_filtered = ch_results_genewise_filtered  // channel: [ tsv ]
    adjacency                 = ch_adjacency                  // channel: [ tsv ]
    versions                  = ch_versions                   // channel: [ versions.yml ]
}

//
// Perform differential analysis
//
include { PROPR_PROPD as PROPD } from "../../../modules/local/propr/propd/main.nf"
include { DESEQ2_DIFFERENTIAL as DESEQ2 } from '../../../modules/nf-core/deseq2/differential/main'
include { FILTER_DIFFTABLE as FILTER_DESEQ2 } from '../../../modules/local/filter_difftable'

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
    ch_control_features
    ch_transcript_lengths

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

    ch_counts
        .join(ch_samplesheet)
        .combine( ch_contrasts )
        .combine( ch_tools_single.deseq2 )
        .unique()
        .multiMap { meta_data, counts, samplesheet, meta_contrast, contrast_variable, reference, target, meta_tools ->
            def meta = meta_data.clone() + ['contrast': meta_contrast.id, 'blocking': meta_contrast.blocking, 'exclude_samples_col': meta_contrast.exclude_samples_col, 'exclude_samples_values': meta_contrast.exclude_samples_values] + meta_tools.clone()
            contrast : [ meta, contrast_variable, reference, target ]
            counts : [ meta, samplesheet, counts ]
        }
        .set { ch_deseq2 }

    // run DESeq2
    DESEQ2 (
        ch_deseq2.contrast,
        ch_deseq2.counts,
        ch_control_features,
        ch_transcript_lengths
    )

    // filter significant DESeq2 results
    // TODO these parameters can be put inside the ext.args of the FILTER_DIFFTABLE module instead,
    // avoiding the explicit channeling of these parameters.
    // To do so, one would need to modify the FILTER_DIFFTABLE module.
    ch_logfc = Channel.value([ params.differential_fc_column, params.differential_min_fold_change ])
    ch_padj = Channel.value([ params.differential_qval_column, params.differential_max_qval ])
    FILTER_DESEQ2(
        DESEQ2.out.results,
        ch_logfc,
        ch_padj
    )

    ch_results_genewise = ch_results_genewise.mix(DESEQ2.out.results)
    ch_results_genewise_filtered = ch_results_genewise_filtered.mix(FILTER_DESEQ2.out.filtered)

    // ----------------------------------------------------
    // Perform differential analysis with LIMMA
    // ----------------------------------------------------

    // TODO: add limma here

    emit:
    results_pairwise          = ch_results_pairwise
    results_pairwise_filtered = ch_results_pairwise_filtered
    results_genewise          = ch_results_genewise
    results_genewise_filtered = ch_results_genewise_filtered
    adjacency                 = ch_adjacency
}

#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/differentialabundance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/differentialabundance
    Website: https://nf-co.re/differentialabundance
    Slack  : https://nfcore.slack.com/channels/differentialabundance
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIALABUNDANCE   } from './workflows/differentialabundance'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.gtf = getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_DIFFERENTIALABUNDANCE {

    take:
    paramsets

    main:

    //
    // WORKFLOW: Run pipeline
    //
    DIFFERENTIALABUNDANCE (
        paramsets
    )

    emit:

    // Preprocessing: Affy
    affy_cel_files             = DIFFERENTIALABUNDANCE.out.affy_cel_files
    affy_raw_expression        = DIFFERENTIALABUNDANCE.out.affy_raw_expression
    affy_norm_expression       = DIFFERENTIALABUNDANCE.out.affy_norm_expression
    affy_annotation            = DIFFERENTIALABUNDANCE.out.affy_annotation
    affy_raw_rds               = DIFFERENTIALABUNDANCE.out.affy_raw_rds

    // Preprocessing: Proteus
    proteus_raw                = DIFFERENTIALABUNDANCE.out.proteus_raw
    proteus_norm               = DIFFERENTIALABUNDANCE.out.proteus_norm
    proteus_plots              = DIFFERENTIALABUNDANCE.out.proteus_plots
    proteus_raw_rdata          = DIFFERENTIALABUNDANCE.out.proteus_raw_rdata
    proteus_norm_rdata         = DIFFERENTIALABUNDANCE.out.proteus_norm_rdata
    proteus_session_info       = DIFFERENTIALABUNDANCE.out.proteus_session_info

    // Preprocessing: GEO
    geo_expression             = DIFFERENTIALABUNDANCE.out.geo_expression
    geo_annotation             = DIFFERENTIALABUNDANCE.out.geo_annotation
    geo_rds                    = DIFFERENTIALABUNDANCE.out.geo_rds

    // Preprocessing: GTF
    gtf_annotation             = DIFFERENTIALABUNDANCE.out.gtf_annotation

    // Differential
    diff_results               = DIFFERENTIALABUNDANCE.out.diff_results
    diff_results_filtered      = DIFFERENTIALABUNDANCE.out.diff_results_filtered
    diff_normalised_matrix     = DIFFERENTIALABUNDANCE.out.diff_normalised_matrix
    diff_variance_stabilised   = DIFFERENTIALABUNDANCE.out.diff_variance_stabilised
    diff_size_factors          = DIFFERENTIALABUNDANCE.out.diff_size_factors
    diff_dispersion_plot       = DIFFERENTIALABUNDANCE.out.diff_dispersion_plot
    diff_md_plot               = DIFFERENTIALABUNDANCE.out.diff_md_plot
    diff_rdata                 = DIFFERENTIALABUNDANCE.out.diff_rdata
    diff_session_info          = DIFFERENTIALABUNDANCE.out.diff_session_info
    diff_annotated             = DIFFERENTIALABUNDANCE.out.diff_annotated

    // Functional: GSEA
    gsea_report_tsv            = DIFFERENTIALABUNDANCE.out.gsea_report_tsv
    gsea_report_html           = DIFFERENTIALABUNDANCE.out.gsea_report_html
    gsea_index_html            = DIFFERENTIALABUNDANCE.out.gsea_index_html
    gsea_heat_map_corr_plot    = DIFFERENTIALABUNDANCE.out.gsea_heat_map_corr_plot
    gsea_ranked_gene_list      = DIFFERENTIALABUNDANCE.out.gsea_ranked_gene_list
    gsea_gene_set_sizes        = DIFFERENTIALABUNDANCE.out.gsea_gene_set_sizes
    gsea_histogram             = DIFFERENTIALABUNDANCE.out.gsea_histogram
    gsea_heatmap               = DIFFERENTIALABUNDANCE.out.gsea_heatmap
    gsea_pvalues_vs_nes_plot   = DIFFERENTIALABUNDANCE.out.gsea_pvalues_vs_nes_plot
    gsea_ranked_list_corr      = DIFFERENTIALABUNDANCE.out.gsea_ranked_list_corr
    gsea_butterfly_plot        = DIFFERENTIALABUNDANCE.out.gsea_butterfly_plot
    gsea_gene_set_tsv          = DIFFERENTIALABUNDANCE.out.gsea_gene_set_tsv
    gsea_gene_set_html         = DIFFERENTIALABUNDANCE.out.gsea_gene_set_html
    gsea_gene_set_heatmap      = DIFFERENTIALABUNDANCE.out.gsea_gene_set_heatmap
    gsea_gene_set_enplot       = DIFFERENTIALABUNDANCE.out.gsea_gene_set_enplot
    gsea_gene_set_dist         = DIFFERENTIALABUNDANCE.out.gsea_gene_set_dist
    gsea_snapshot              = DIFFERENTIALABUNDANCE.out.gsea_snapshot
    gsea_archive               = DIFFERENTIALABUNDANCE.out.gsea_archive
    gsea_rpt                   = DIFFERENTIALABUNDANCE.out.gsea_rpt

    // Functional: gprofiler2
    gprofiler2_html            = DIFFERENTIALABUNDANCE.out.gprofiler2_html
    gprofiler2_all_enrichment  = DIFFERENTIALABUNDANCE.out.gprofiler2_all_enrichment
    gprofiler2_sub_enrichment  = DIFFERENTIALABUNDANCE.out.gprofiler2_sub_enrichment
    gprofiler2_plot_png        = DIFFERENTIALABUNDANCE.out.gprofiler2_plot_png
    gprofiler2_sub_plot        = DIFFERENTIALABUNDANCE.out.gprofiler2_sub_plot
    gprofiler2_rds             = DIFFERENTIALABUNDANCE.out.gprofiler2_rds
    gprofiler2_filtered_gmt    = DIFFERENTIALABUNDANCE.out.gprofiler2_filtered_gmt

    // Functional: decoupler
    decoupler_estimate         = DIFFERENTIALABUNDANCE.out.decoupler_estimate
    decoupler_pvals            = DIFFERENTIALABUNDANCE.out.decoupler_pvals
    decoupler_png              = DIFFERENTIALABUNDANCE.out.decoupler_png

    // Functional: common
    functional_session_info    = DIFFERENTIALABUNDANCE.out.functional_session_info

    // Plotting
    plot_exploratory           = DIFFERENTIALABUNDANCE.out.plot_exploratory
    plot_volcanos              = DIFFERENTIALABUNDANCE.out.plot_volcanos

    // ShinyNGS
    shinyngs_data              = DIFFERENTIALABUNDANCE.out.shinyngs_data
    shinyngs_app_file          = DIFFERENTIALABUNDANCE.out.shinyngs_app_file

    // Report
    report_html                = DIFFERENTIALABUNDANCE.out.report_html
    report_bundle              = DIFFERENTIALABUNDANCE.out.report_bundle

    // Versions
    nfcore_versions            = DIFFERENTIALABUNDANCE.out.nfcore_versions
    collated_versions          = DIFFERENTIALABUNDANCE.out.collated_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //

    NFCORE_DIFFERENTIALABUNDANCE (
        PIPELINE_INITIALISATION.out.paramsets
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
    )

    //
    // Build category channels for publishing (name-tagging happens here)
    //
    def out = NFCORE_DIFFERENTIALABUNDANCE.out

    ch_pub_preprocessing = out.affy_cel_files.map { it -> ['affy_cel_files', it[0], it[1..-1]] }
        .mix(out.affy_raw_expression.map        { meta, file -> ['affy_raw_expression', meta, file] })
        .mix(out.affy_norm_expression.map       { meta, file -> ['affy_norm_expression', meta, file] })
        .mix(out.affy_annotation.map            { meta, file -> ['affy_annotation', meta, file] })
        .mix(out.affy_raw_rds.map               { meta, file -> ['affy_raw_rds', meta, file] })
        .mix(out.proteus_raw.map                { meta, file -> ['proteus_raw', meta, file] })
        .mix(out.proteus_norm.map               { meta, file -> ['proteus_norm', meta, file] })
        .mix(out.proteus_plots.map              { meta, file -> ['proteus_plots', meta, file] })
        .mix(out.proteus_raw_rdata.map          { meta, file -> ['proteus_raw_rdata', meta, file] })
        .mix(out.proteus_norm_rdata.map         { meta, file -> ['proteus_norm_rdata', meta, file] })
        .mix(out.proteus_session_info.map       { meta, file -> ['proteus_session_info', meta, file] })
        .mix(out.geo_expression.map             { meta, file -> ['geo_expression', meta, file] })
        .mix(out.geo_annotation.map             { meta, file -> ['geo_annotation', meta, file] })
        .mix(out.geo_rds.map                    { meta, file -> ['geo_rds', meta, file] })
        .mix(out.gtf_annotation.map             { meta, file -> ['gtf_annotation', meta, file] })

    ch_pub_differential = out.diff_results.map              { _key, meta, file -> ['results', meta, file] }
        .mix(out.diff_results_filtered.map  { _key, meta, file -> ['results_filtered', meta, file] })
        .mix(out.diff_normalised_matrix.map { meta, file -> ['normalised_matrix', meta, file] })
        .mix(out.diff_variance_stabilised.map { meta, file -> ['variance_stabilised_matrix', meta, file] })
        .mix(out.diff_size_factors.map      { meta, file -> ['size_factors', meta, file] })
        .mix(out.diff_dispersion_plot.map   { meta, file -> ['dispersion_plot', meta, file] })
        .mix(out.diff_md_plot.map           { meta, file -> ['md_plot', meta, file] })
        .mix(out.diff_rdata.map             { meta, file -> ['rdata', meta, file] })
        .mix(out.diff_session_info.map      { meta, file -> ['session_info', meta, file] })
        .mix(out.diff_annotated.map         { meta, file -> ['annotated', meta, file] })

    ch_pub_functional = out.gsea_report_tsv.map          { meta, ref, target -> ['gsea_report_tsv', meta, [ref, target]] }
        .mix(out.gsea_report_html.map       { meta, ref, target -> ['gsea_report_html', meta, [ref, target]] })
        .mix(out.gsea_index_html.map        { meta, file -> ['gsea_index_html', meta, file] })
        .mix(out.gsea_heat_map_corr_plot.map { meta, file -> ['gsea_heat_map_corr_plot', meta, file] })
        .mix(out.gsea_ranked_gene_list.map  { meta, file -> ['gsea_ranked_gene_list', meta, file] })
        .mix(out.gsea_gene_set_sizes.map    { meta, file -> ['gsea_gene_set_sizes', meta, file] })
        .mix(out.gsea_histogram.map         { meta, file -> ['gsea_histogram', meta, file] })
        .mix(out.gsea_heatmap.map           { meta, file -> ['gsea_heatmap', meta, file] })
        .mix(out.gsea_pvalues_vs_nes_plot.map { meta, file -> ['gsea_pvalues_vs_nes_plot', meta, file] })
        .mix(out.gsea_ranked_list_corr.map  { meta, file -> ['gsea_ranked_list_corr', meta, file] })
        .mix(out.gsea_butterfly_plot.map    { meta, file -> ['gsea_butterfly_plot', meta, file] })
        .mix(out.gsea_gene_set_tsv.map      { meta, file -> ['gsea_gene_set_tsv', meta, file] })
        .mix(out.gsea_gene_set_html.map     { meta, file -> ['gsea_gene_set_html', meta, file] })
        .mix(out.gsea_gene_set_heatmap.map  { meta, file -> ['gsea_gene_set_heatmap', meta, file] })
        .mix(out.gsea_gene_set_enplot.map   { meta, file -> ['gsea_gene_set_enplot', meta, file] })
        .mix(out.gsea_gene_set_dist.map     { meta, file -> ['gsea_gene_set_dist', meta, file] })
        .mix(out.gsea_snapshot.map          { meta, file -> ['gsea_snapshot', meta, file] })
        .mix(out.gsea_archive.map           { meta, file -> ['gsea_archive', meta, file] })
        .mix(out.gsea_rpt.map               { meta, file -> ['gsea_rpt', meta, file] })
        .mix(out.gprofiler2_html.map        { meta, file -> ['gprofiler2_html', meta, file] })
        .mix(out.gprofiler2_all_enrichment.map { meta, file -> ['gprofiler2_all_enrichment', meta, file] })
        .mix(out.gprofiler2_sub_enrichment.map { meta, file -> ['gprofiler2_sub_enrichment', meta, file] })
        .mix(out.gprofiler2_plot_png.map    { meta, file -> ['gprofiler2_plot_png', meta, file] })
        .mix(out.gprofiler2_sub_plot.map    { meta, file -> ['gprofiler2_sub_plot', meta, file] })
        .mix(out.gprofiler2_rds.map         { meta, file -> ['gprofiler2_rds', meta, file] })
        .mix(out.gprofiler2_filtered_gmt.map { meta, file -> ['gprofiler2_filtered_gmt', meta, file] })
        .mix(out.decoupler_estimate.map     { meta, file -> ['decoupler_estimate', meta, file] })
        .mix(out.decoupler_pvals.map        { meta, file -> ['decoupler_pvals', meta, file] })
        .mix(out.decoupler_png.map          { meta, file -> ['decoupler_png', meta, file] })
        .mix(out.functional_session_info.map { meta, file -> ['session_info', meta, file] })

    ch_pub_plotting = out.plot_exploratory.map { meta, file -> ['exploratory', meta, file] }
        .mix(out.plot_volcanos.map          { meta, file -> ['differential_volcanos', meta, file] })

    ch_pub_shinyngs = out.shinyngs_data.map { meta, file -> ['shinyngs_data', meta, file] }
        .mix(out.shinyngs_app_file.map      { meta, file -> ['shinyngs_app', meta, file] })

    ch_pub_report = out.report_html.map { meta, file -> ['report_html', meta, file] }
        .mix(out.report_bundle.map          { meta, file -> ['report_bundle', meta, file] })

    ch_pub_versions = out.nfcore_versions.map { file -> ['versions', [:], file] }
        .mix(out.collated_versions.map      { file -> ['collated_versions', [:], file] })

    publish:
    preprocessing = ch_pub_preprocessing
    differential  = ch_pub_differential
    functional    = ch_pub_functional
    plotting      = ch_pub_plotting
    shinyngs      = ch_pub_shinyngs
    report        = ch_pub_report
    versions      = ch_pub_versions
}

output {
    preprocessing {
        path { name, meta, file ->
            def folder = [
                // AFFY
                affy_raw_expression        : 'tables/processed_abundance',
                affy_norm_expression       : 'tables/processed_abundance',
                affy_annotation            : 'tables/annotation',
                affy_raw_rds               : 'other/affy',
                affy_cel_files             : 'untar',

                // PROTEUS
                proteus_raw                : 'tables/proteus',
                proteus_norm               : 'tables/proteus',
                proteus_plots              : 'plots/proteus',
                proteus_raw_rdata          : 'other/proteus',
                proteus_norm_rdata         : 'other/proteus',
                proteus_session_info       : 'other/proteus',

                // GEO SOFT
                geo_expression             : 'tables/processed_abundance',
                geo_annotation             : 'tables/annotation',
                geo_rds                    : 'other/affy',

                // GTF
                gtf_annotation             : 'tables/annotation',

            ][name] ?: name
            def target = (name in ['proteus_plots', 'proteus_raw_rdata', 'proteus_norm_rdata']) \
                ? "${folder}/${meta.paramset_name}/${meta.contrast}/" \
                : "${folder}/${meta.paramset_name}/"
            return target
        }
    }
    differential {
        path { name, meta, file ->
            def folder = [
                results                    : 'tables/differential',
                results_filtered           : 'tables/differential',
                annotated                  : 'tables/differential',
                normalised_matrix          : 'tables/processed_abundance',
                variance_stabilised_matrix : 'tables/processed_abundance',
                size_factors               : "other/${meta.params.differential_method}",
                dispersion_plot            : 'plots/qc',
                md_plot                    : 'plots/qc',
                rdata                      : "other/${meta.params.differential_method}",
                session_info               : "other/${meta.params.differential_method}"
            ][name] ?: name
            file >> "${folder}/${meta.paramset_name}/"
        }
    }
    functional {
        path { name, meta, file ->
            def method = meta.params.functional_method
            def folder = [
                // GSEA
                gsea_report_tsv           : 'report/gsea',
                gsea_report_html          : 'report/gsea',
                gsea_index_html           : 'report/gsea',
                gsea_heat_map_corr_plot   : 'report/gsea',
                gsea_ranked_gene_list     : 'report/gsea',
                gsea_gene_set_sizes       : 'report/gsea',
                gsea_histogram            : 'report/gsea',
                gsea_heatmap              : 'report/gsea',
                gsea_pvalues_vs_nes_plot  : 'report/gsea',
                gsea_ranked_list_corr     : 'report/gsea',
                gsea_butterfly_plot       : 'report/gsea',
                gsea_gene_set_tsv         : 'report/gsea',
                gsea_gene_set_html        : 'report/gsea',
                gsea_gene_set_heatmap     : 'report/gsea',
                gsea_gene_set_enplot      : 'report/gsea',
                gsea_gene_set_dist        : 'report/gsea',
                gsea_snapshot             : 'report/gsea',
                gsea_archive              : 'report/gsea',
                gsea_rpt                  : 'report/gsea',

                // GPROFILER2
                gprofiler2_all_enrichment : 'tables/gprofiler2',
                gprofiler2_sub_enrichment : 'tables/gprofiler2',
                gprofiler2_html           : 'plots/gprofiler2',
                gprofiler2_plot_png       : 'plots/gprofiler2',
                gprofiler2_sub_plot       : 'plots/gprofiler2',
                gprofiler2_rds            : 'other/gprofiler2',
                gprofiler2_filtered_gmt   : 'other/gprofiler2',

                // DECOUPLER
                decoupler_estimate        : 'tables/decoupler',
                decoupler_pvals           : 'tables/decoupler',
                decoupler_png             : 'plots/decoupler',

                // common outputs
                session_info               : "other/${method}"
            ][name] ?: name
            def gene_set_name = (method == 'gsea' && meta.params.gene_sets_files) \
                ? meta.params.gene_sets_files.tokenize('/')[-1].replaceFirst(/\.[^.]+$/, '') \
                : null
            def target = (method == 'gsea') ? "${folder}/${meta.paramset_name}/${meta.id}/${gene_set_name}" \
                : (method == 'gprofiler2') ? "${folder}/${meta.paramset_name}/${meta.id}/" \
                : "${folder}/${meta.paramset_name}/"
            return target
        }
    }
    plotting {
        path { name, meta, file ->
            def folder = [
                exploratory           : 'plots/exploratory',
                differential_volcanos : 'plots/differential',
            ][name] ?: name
            file >> "${folder}/${meta.paramset_name}/"
        }
    }
    shinyngs {
        path { name, meta, file ->
            file >> "shinyngs_app/${meta.paramset_name}/"
        }
    }
    report {
        path { name, meta, file ->
            file >> "report/${meta.paramset_name}/"
        }
    }
    versions {
        path { name, meta, file ->
            file >> "pipeline_info/"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

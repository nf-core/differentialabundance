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
    preprocessing = DIFFERENTIALABUNDANCE.out.preprocessing
    differential  = DIFFERENTIALABUNDANCE.out.differential
    functional    = DIFFERENTIALABUNDANCE.out.functional
    plotting      = DIFFERENTIALABUNDANCE.out.plotting
    shinyngs      = DIFFERENTIALABUNDANCE.out.shinyngs
    report        = DIFFERENTIALABUNDANCE.out.report
    versions      = DIFFERENTIALABUNDANCE.out.versions
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
        params.hook_url,
    )

    publish:
    preprocessing = NFCORE_DIFFERENTIALABUNDANCE.out.preprocessing
    differential  = NFCORE_DIFFERENTIALABUNDANCE.out.differential
    functional    = NFCORE_DIFFERENTIALABUNDANCE.out.functional
    plotting      = NFCORE_DIFFERENTIALABUNDANCE.out.plotting
    shinyngs      = NFCORE_DIFFERENTIALABUNDANCE.out.shinyngs
    report        = NFCORE_DIFFERENTIALABUNDANCE.out.report
    versions      = NFCORE_DIFFERENTIALABUNDANCE.out.versions
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

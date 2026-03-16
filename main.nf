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
    tables        = DIFFERENTIALABUNDANCE.out.tables
    plots         = DIFFERENTIALABUNDANCE.out.plots
    report        = DIFFERENTIALABUNDANCE.out.report
    shinyngs_app  = DIFFERENTIALABUNDANCE.out.shinyngs_app
    other         = DIFFERENTIALABUNDANCE.out.other
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
    tables        = NFCORE_DIFFERENTIALABUNDANCE.out.tables
    plots         = NFCORE_DIFFERENTIALABUNDANCE.out.plots
    report        = NFCORE_DIFFERENTIALABUNDANCE.out.report
    shinyngs_app  = NFCORE_DIFFERENTIALABUNDANCE.out.shinyngs_app
    other         = NFCORE_DIFFERENTIALABUNDANCE.out.other
    versions      = NFCORE_DIFFERENTIALABUNDANCE.out.versions
}

// TODO organize output folders to meta.paramset_name/tables/folder instead later on,
// but for now we do this so that snapshots are comparable.
output {
    tables {
        path { name, meta, file ->
            def folder = [
                affy_raw_expression    : 'processed_abundance',
                affy_norm_expression   : 'processed_abundance',
                affy_annotation        : 'annotation',
                proteus_raw            : 'proteus',
                proteus_norm           : 'proteus',
                geo_expression         : 'processed_abundance',
                geo_annotation         : 'annotation',
                gtf_annotation         : 'annotation',
                immunedeconv           : 'immunedeconv',
                differential_results          : 'differential',
                differential_results_filtered : 'differential',
                normalised_matrix             : 'processed_abundance',
                variance_stabilised_matrix    : 'processed_abundance',
                differential_annotated        : 'differential',
                gprofiler2_all_enrichment     : 'gprofiler2',
                gprofiler2_sub_enrichment     : 'gprofiler2',
                decoupler_estimate            : 'decoupler',
                decoupler_pvals               : 'decoupler',
            ][name] ?: name
            file >> "tables/${folder}/${meta.paramset_name}/"
        }
    }
    plots {
        path { name, meta, file ->
            def folder = [
                proteus_dendrogram              : 'proteus',
                proteus_mean_variance           : 'proteus',
                proteus_raw_distributions       : 'proteus',
                proteus_normalized_distributions: 'proteus',
                immunedeconv                    : 'immunedeconv',
                differential_qc                 : 'qc',
                gprofiler2_plot                 : 'gprofiler2',
                gprofiler2_sub_plot             : 'gprofiler2',
                decoupler                       : 'decoupler',
                exploratory                     : 'exploratory',
                differential_volcanos           : 'differential',
            ][name] ?: name
            file >> "plots/${folder}/${meta.paramset_name}/"
        }
    }
    shinyngs_app {
        path { name, meta, file ->
            file >> "shinyngs_app/${meta.paramset_name}/"
        }
    }
    other {
        path { name, meta, file ->
            def folder = [
                affy_raw_rds          : 'affy',
                proteus_raw_rdata     : 'proteus',
                proteus_norm_rdata    : 'proteus',
                proteus_session_info  : 'proteus',
                differential_other    : meta.method_differential,
                gprofiler2_other      : 'gprofiler2',
                geo_rds               : 'geo',
            ][name] ?: name
            file >> "other/${folder}/${meta.paramset_name}/"
        }
    }
    report {
        path { name, meta, file ->
            def folder = [
                gprofiler2    : 'gprofiler2',
                gsea          : 'gsea',
                report_html   : '',
                report_bundle : '',
            ][name] ?: name
            def subpath = folder ? "report/${folder}/${meta.paramset_name}/" : "report/${meta.paramset_name}/"
            file >> subpath
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

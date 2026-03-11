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

    DIFFERENTIALABUNDANCE (
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
    tables = DIFFERENTIALABUNDANCE.out.tables
    plots = DIFFERENTIALABUNDANCE.out.plots
    report = DIFFERENTIALABUNDANCE.out.report
    shinyngs_app = DIFFERENTIALABUNDANCE.out.shinyngs_app
    other = DIFFERENTIALABUNDANCE.out.other
    pipeline_info = DIFFERENTIALABUNDANCE.out.pipeline_info
}

output {
    tables {
        path { subdir, file -> "tables/${subdir}/${file.name}" }
    }

    plots {
        path { subdir, file -> "plots/${subdir}/${file.name}" }
    }

    report {
        path { subdir, file -> "report/${subdir}/${file.name}" }
    }

    shinyngs_app {
        path { subdir, file -> "shinyngs_app/${subdir}/${file.name}" }
    }

    other {
        path { subdir, file -> "other/${subdir}/${file.name}" }
    }

    pipeline_info {
        path { subdir, file -> "${subdir}/${file.name}" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDifferentialabundance.initialise(params, log)

def checkPathParamList = [ params.input, params.gtf, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../modules/nf-core/shinyngs/staticexploratory/main'  
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../modules/nf-core/shinyngs/staticdifferential/main'  
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'  
include { DESEQ2_DIFFERENTIAL                               } from '../modules/nf-core/deseq2/differential/main'    
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DIFFERENTIALABUNDANCE {

    if ( params.gtf.endsWith('.gz') ){
        file_gtf_in = file(params.gtf)
        GUNZIP_GTF([["id": file_gtf_in.simpleName], file_gtf_in])
        file_gtf = GUNZIP_GTF.out.gunzip
    } else{
        file_gtf = file(params.gtf) 
    }

    exp_meta = [ "id": params.study_name  ]

    // Get a features table from the GTF and combine with the matrix and sample
    // annotation (fom = features/ observations/ matrix)

    GTF_TO_TABLE( file_gtf, [[ "id":""], []])
    
    // Combine features with the observations and matrices to create a FOM
    // where these things can travel together

    ch_fom = GTF_TO_TABLE.out.feature_annotation
        .map{
            tuple( exp_meta, file(params.input), it[1], file(params.matrix))
        }
   
    // Channel for the contrasts file
    
    ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])

    // Check compatibility of FOM elements and contrasts

    VALIDATOR(  
        ch_fom,
        ch_contrasts_file
    )
 
    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately

    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )

    // Run the DESeq differential module, which doesn't take the feature
    // annotations 

    ch_samples_and_matrix = VALIDATOR.out.fom.map{
        tuple(it[1], it[3])
    }

    DESEQ2_DIFFERENTIAL (
        ch_contrasts.combine(ch_samples_and_matrix)
    )

    // Let's make the simplifying assumption that the processed matrices from
    // the DESeq runs are the same across contrasts. We run the DESeq process
    // with matrices once for each contrast because DESeqDataSetFromMatrix()
    // takes the model, and the model can vary between contrasts if the
    // blocking factors included differ. But the normalised and
    // variance-stabilised matrices are not (IIUC) impacted by the model.

    ch_processed_matrices = DESEQ2_DIFFERENTIAL.out.normalised_counts
        .join(DESEQ2_DIFFERENTIAL.out.vst_counts)
        .map{ it.tail() }
        .first()
   
    // The exploratory plots are made by coloring by every unique variable used
    // to define contrasts

    ch_contrast_variables = ch_contrasts
        .map{
            [ "id": it.variable ]
        }
        .unique()

    ch_fom_plot_inputs = VALIDATOR.out.fom
        .combine(ch_processed_matrices)                         // Add processed marices to what we have in the FOM
        .map{
            tuple(it[0], it[1], it[2], [ it[3], it[4], it[5] ]) // Remove the experiment meta and group the matrices
        }
 
    PLOT_EXPLORATORY(
        ch_contrast_variables
            .combine(ch_fom_plot_inputs.map{ it.tail() })
    )

    // Differential analysis using the results of DESeq2

    PLOT_DIFFERENTIAL(
        DESEQ2_DIFFERENTIAL.out.results, 
        ch_fom_plot_inputs.first()
    )

    // Gather software versions

    ch_versions = GTF_TO_TABLE.out.versions
        .mix(VALIDATOR.out.versions)
        .mix(DESEQ2_DIFFERENTIAL.out.versions)
        .mix(PLOT_EXPLORATORY.out.versions)
        .mix(PLOT_DIFFERENTIAL.out.versions)

    if ( params.gtf.endsWith('.gz') ){
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
//    workflow_summary    = WorkflowDifferentialabundance.paramsSummaryMultiqc(workflow, summary_params)
//    ch_workflow_summary = Channel.value(workflow_summary)

//    methods_description    = WorkflowDifferentialabundance.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
//    ch_methods_description = Channel.value(methods_description)

//    ch_multiqc_files = Channel.empty()
//    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
//    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

//    MULTIQC (
 //       ch_multiqc_files.collect(),
//        ch_multiqc_config.collect().ifEmpty([]),
//        ch_multiqc_custom_config.collect().ifEmpty([]),
//        ch_multiqc_logo.collect().ifEmpty([])
//    )
//    multiqc_report = MULTIQC.out.report.toList()
//    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

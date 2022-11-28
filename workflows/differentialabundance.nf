/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDifferentialabundance.initialise(params, log)

def checkPathParamList = [ params.input, params.gtf ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.report_file) { ch_report_file = file(params.report_file) } else { ch_report_file = [] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXECUTE_REPORT                                    } from '../modules/local/executereport/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
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

    ch_versions = Channel.empty()
   
    // Get feature annotations from a GTF file, gunzip if necessary
 
    file_gtf_in = file(params.gtf)
    file_gtf = [ [ "id": file_gtf_in.simpleName ], file_gtf_in ] 

    if ( params.gtf.endsWith('.gz') ){
        GUNZIP_GTF(file_gtf)
        file_gtf = GUNZIP_GTF.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
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
    // downstream plots separately. 
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )
        .map{
            it.blocking = it.blocking.replace('NA', '')
            it
        }

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

    // Output report

    if (params.report_file) {
        EXECUTE_REPORT(
            ch_report_file,
            params
        )
    }

    // Gather software versions

    ch_versions = ch_versions
        .mix(GTF_TO_TABLE.out.versions)
        .mix(VALIDATOR.out.versions)
        .mix(DESEQ2_DIFFERENTIAL.out.versions)
        .mix(PLOT_EXPLORATORY.out.versions)
        .mix(PLOT_DIFFERENTIAL.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Validate input parameters
// TODO make it specific for logratioanalysis || or maybe it common for the whole workflow
WorkflowDifferentialabundance.initialise(params,log)

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input,
    params.matrix,
    params.gtf
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Define meta
def exp_meta = [ "id": params.study_name  ]

// Check input samplesheet
if (params.input) { 
    input_file = file(params.input, checkIfExists: true)
    ch_input = Channel.of([ exp_meta, input_file ]) 
    } else { 
        exit 1, 'Input samplesheet not specified!' 
    }
// Check input matrix
if (params.matrix) { 
        matrix_file = file(params.matrix, checkIfExists: true)
        ch_in_raw = Channel.of([ exp_meta, matrix_file])
    } else { 
        error("Input matrix not specified!")
    }

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

include { PROPR_PROPR                                       } from '../modules/local/propr/propr/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { PROPR_LOGRATIO                                    } from '../modules/nf-core/propr/logratio/main'
include { PROPR_PROPD                                       } from '../modules/nf-core/propr/propd/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow LOGRATIOANALYSIS {


    /*
     * DATA HANDLING
     * Fetch and parse the dataset and metadata.
     * Filter data if needed.
     * Replace zeros if needed.
     */

    // Set up some basic variables
    ch_versions = Channel.empty()

    // matrix
    ch_matrix = ch_in_raw

    // Firstly Filter the input matrix
    // TODO filter both genes and samples(cells) according to some criteria
    CUSTOM_MATRIXFILTER(
        ch_matrix,
        ch_input
    )

    // TODO currently the zeros can be directly handled in the downstream analysis by replacing them with the min value
    // TODO add zero imputation methods

    // TODO need to handle different data structures here (eg. affy, etc)


    /*
     * LOGRATIO-BASED CORRELATION ANALYSIS
     * Compute basis shrinkage partial correlation, or proportionality
     */
    PROPR_LOGRATIO(
        CUSTOM_MATRIXFILTER.out.filtered
    )


    /*
    * DIFFERENTIAL PROPORTIONALITY ANALYSIS
    * Compute differential proportionality coefficients
    */
    PROPR_PROPD(
        CUSTOM_MATRIXFILTER.out.filtered
    )

    // TODO add GREA


}
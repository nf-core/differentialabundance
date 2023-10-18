/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowDifferentialabundance.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
def exp_meta = [ "id": params.study_name  ]
if (params.input) { 
    ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ]) 
} else { 
    exit 1, 'Input samplesheet not specified!' 
}

// handle different data formats
if (params.study_type == 'affy_array'){

    if (params.affy_cel_files_archive) {
        ch_celfiles = Channel.of([ exp_meta, file(params.affy_cel_files_archive, checkIfExists: true) ])
    } else {
        error("CEL files archive not specified!")
    }
    
} else {
    // If this is not microarray data or maxquant output, and this an RNA-seq dataset,
    // then assume we're reading from a matrix

    if (params.study_type == "rnaseq" && params.matrix) {
        matrix_file = file(params.matrix, checkIfExists: true)
        ch_in_raw = Channel.of([ exp_meta, matrix_file])
    } else {
        error("Input matrix not specified!")
    }

}

report_file    = file(params.report_file, checkIfExists: true)
logo_file      = file(params.logo_file, checkIfExists: true)
css_file       = file(params.css_file, checkIfExists: true)
citations_file = file(params.citations_file, checkIfExists: true)


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

include { PROPR_PROPR as PROPR_PARTIALCORRELATION           } from '../modules/local/propr/propr/main'
include { PROPR_PROPR as PROPR_PROPORTIONALITY              } from '../modules/local/propr/propr/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR                                             } from '../modules/nf-core/untar/main.nf'
include { AFFY_JUSTRMA as AFFY_JUSTRMA_RAW                  } from '../modules/nf-core/affy/justrma/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { PROPR_LOGRATIO                                    } from '../modules/nf-core/propr/logratio/main'
include { PROPR_PROPD as PROPR_DIFFERENTIALPROPORTIONALITY  } from '../modules/nf-core/propr/propd/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOGRATIOANALYSIS {

    // Set up some basic variables
    ch_versions = Channel.empty()

    /*
     * GET DATA
     * it parses the expression matrix, the features metadata, and the contrast file
     */

    // format data
    if (params.study_type == 'affy_array'){

        // Uncompress the CEL files archive
        UNTAR ( ch_celfiles )
        ch_affy_input = ch_input
            .join(UNTAR.out.untar)

        // Run affy to derive the matrix
        AFFY_JUSTRMA_RAW (
            ch_affy_input,
            [[],[]]
        )

        // Fetch affy outputs and reset the meta
        ch_in_raw = AFFY_JUSTRMA_RAW.out.expression
        ch_affy_platform_features = AFFY_JUSTRMA_RAW.out.annotation
        
        ch_versions = ch_versions
            .mix(UNTAR.out.versions)
            .mix(AFFY_JUSTRMA_RAW.out.versions)
    }

    // If user has provided a feature annotation table, use that
    if (params.features){
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)])

    } else if (params.study_type == 'affy_array'){
        ch_features = ch_affy_platform_features

    } else if (params.study_type == 'rnaseq'){

        // Otherwise we can just use the matrix input; save it to the workdir so that it does not
        // just appear wherever the user runs the pipeline
        matrix_as_anno_filename = "${workflow.workDir}/matrix_as_anno.${matrix_file.getExtension()}"
        ch_features = ch_in_raw
            .map{ meta, matrix ->
                matrix.copyTo(matrix_as_anno_filename)
                [ meta, file(matrix_as_anno_filename) ]
            }
    }

    // Channel for the contrasts file
    ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])

    /* 
     * VALIDATE DATA
     * It checks the matrix, the metadata and contrast files are consistent between them.
     * It also parses the matrix, so that the rows are named by the observations, 
     * and the columns are named by the features. The rest of metadata are kept separately.
     */

    VALIDATOR(
        ch_input.join(ch_in_raw),
        ch_features,
        ch_contrasts_file
    ) 
    ch_matrix = VALIDATOR.out.assays
    ch_versions = ch_versions
            .mix(VALIDATOR.out.versions)

    // TODO use contrast

    // TODO logratio analysis should be applied on raw read counts. 
    // We may want to do some input data checking related on this

    /* 
     * FILTER DATA
     * filter lowly expressed genes
     */

    CUSTOM_MATRIXFILTER(
        ch_matrix,
        ch_input
    )
    ch_matrix = CUSTOM_MATRIXFILTER.out.filtered
    ch_versions = ch_versions
            .mix(CUSTOM_MATRIXFILTER.out.versions)

    /*
     * PROCESS DATA
     * Replace zeros if needed.
     * Normalize data if needed.
     */

    // TODO currently the zeros can be directly handled in the downstream analysis by replacing them with the min value
    // TODO add zero imputation methods

    // TODO add normalization blocks

    /*
     * LOGRATIO-BASED CORRELATION ANALYSIS
     * Compute basis shrinkage partial correlation, or proportionality
     */

    PROPR_PARTIALCORRELATION(
        ch_matrix
    )
    PROPR_PROPORTIONALITY(
        ch_matrix
    )
    ch_versions = ch_versions
            .mix(PROPR_PARTIALCORRELATION.out.versions)
            .mix(PROPR_PROPORTIONALITY.out.versions)

    /*
    * DIFFERENTIAL PROPORTIONALITY ANALYSIS
    * Compute differential proportionality coefficients
    */
    PROPR_DIFFERENTIALPROPORTIONALITY(
        ch_matrix,
        ch_input
    )
    ch_versions = ch_versions
            .mix(PROPR_DIFFERENTIALPROPORTIONALITY.out.versions)

    // TODO add GREA

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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

include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { PROPR_LOGRATIO                                    } from '../modules/nf-core/propr/logratio/main'
include { PROPR_PROPD as PROPR_DIFFERENTIALPROPORTIONALITY  } from '../modules/nf-core/propr/propd/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOGRATIOANALYSIS {

    take:
    ch_input
    ch_in_raw


    main:

    /*
     * DATA HANDLING
     * Fetch and parse the dataset and metadata.
     * Filter data if needed.
     * Replace zeros if needed.
     * Normalize data if needed.
     */

    // Set up some basic variables
    ch_versions = Channel.empty()

    // matrix
    ch_matrix = ch_in_raw

    // Firstly Filter the input matrix
    CUSTOM_MATRIXFILTER(
        ch_matrix,
        ch_input
    )

    // TODO need to handle different data structures here (eg. affy, etc)

    // TODO use contrast

    // TODO currently the zeros can be directly handled in the downstream analysis by replacing them with the min value
    // TODO add zero imputation methods

    // TODO add normalization blocks


    /*
     * LOGRATIO-BASED CORRELATION ANALYSIS
     * Compute basis shrinkage partial correlation, or proportionality
     */

    // compute partial correlation with shrinkage
    PROPR_PARTIALCORRELATION(
        CUSTOM_MATRIXFILTER.out.filtered
    )

    // compute proportionality
    PROPR_PROPORTIONALITY(
        CUSTOM_MATRIXFILTER.out.filtered
    )


    /*
    * DIFFERENTIAL PROPORTIONALITY ANALYSIS
    * Compute differential proportionality coefficients
    */
    PROPR_DIFFERENTIALPROPORTIONALITY(
        CUSTOM_MATRIXFILTER.out.filtered,
        ch_input
    )

    // TODO add GREA


}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set mandatory input channel
def exp_meta = [ "id": params.study_name  ]
ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ])

// Handle different data formats
// TODO now it only handles affy_array and rnaseq. It should be able to handle all the data types that differentialabundance.nf handles.
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

// TODO need to make the containers and put these modules into nf-core
include { PROPR_PROPR as PROPR_PARTIALCORRELATION           } from '../modules/local/propr/propr/main'
include { PROPR_PROPR as PROPR_PROPORTIONALITY              } from '../modules/local/propr/propr/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                                             } from '../modules/nf-core/untar/main.nf'
include { AFFY_JUSTRMA                                      } from '../modules/nf-core/affy/justrma/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { PROPR_PROPD as PROPR_DIFFERENTIALPROPORTIONALITY  } from '../modules/nf-core/propr/propd/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
// TODO need to implement the code to record the metrics into multiqc report
def multiqc_report = []

workflow LOGRATIOANALYSIS {

    // Set up some basic variables
    ch_versions = Channel.empty()

    /*
     * GET DATA
     * it parses the expression matrix, the features metadata, and the contrast file
     */

    // Get expression matrix
    // This is required when data type is not rnaseq
    // TODO add code for other data types
    if (params.study_type == 'affy_array'){

        // Uncompress the CEL files archive

        UNTAR ( ch_celfiles )

        ch_affy_input = ch_input
            .join(UNTAR.out.untar)

        // Run affy to derive the matrix

        AFFY_JUSTRMA (
            ch_affy_input,
            [[],[]]
        )

        // Fetch affy outputs and reset the meta

        ch_in_raw = AFFY_JUSTRMA.out.expression
        ch_affy_platform_features = AFFY_JUSTRMA.out.annotation
        
        ch_versions = ch_versions
            .mix(UNTAR.out.versions)
            .mix(AFFY_JUSTRMA.out.versions)
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
        VALIDATOR.out.sample_meta
    )

    ch_input = VALIDATOR.out.sample_meta
    ch_matrix = CUSTOM_MATRIXFILTER.out.filtered
    ch_samples_and_matrix = ch_input.join(ch_matrix)     // -> meta, samplesheet, filtered matrix

    ch_versions = ch_versions
            .mix(CUSTOM_MATRIXFILTER.out.versions)

    // TODO run and check what is VALIDATOR.out.sample_meta

    /*
     * HANDLE ZEROS
     */

    // TODO currently the zeros can be directly handled in the downstream analysis by replacing them with the min value
    // TODO add zero imputation methods

    /*
     * LOGRATIO-BASED ANALYSIS
     */

    // TODO should I also put the samplesheet information into logratio modules?

    // compute logratio partial correlation with basis shrinkage
    if (params.run_partial_correlation){
        PROPR_PARTIALCORRELATION(
            ch_matrix
        )
        ch_versions = ch_versions
            .mix(PROPR_PARTIALCORRELATION.out.versions)
    }

    // compute proportionality
    if (params.run_proportionality){
        PROPR_PROPORTIONALITY(
            ch_matrix
        )
        ch_versions = ch_versions
                .mix(PROPR_PROPORTIONALITY.out.versions)
    }

    // compute differential proportionality
    if (params.run_differential_proportionality){
        PROPR_DIFFERENTIALPROPORTIONALITY(
            ch_matrix,
            ch_input
        )
        ch_versions = ch_versions
                .mix(PROPR_DIFFERENTIALPROPORTIONALITY.out.versions)

        // TODO add GREA
    }

    // TODO plot exploratory figures
    // TODO add results analysis, etc

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // TODO handle report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

// Check optional parameters
if (params.control_features) { ch_control_features = file(params.control_features, checkIfExists: true) } else { ch_control_features = [[],[]] } 
if (params.gsea_run) { gene_sets_file = file(params.gsea_gene_sets, checkIfExists: true) } else { gene_sets_file = [] } 

report_file = file(params.report_file, checkIfExists: true)
logo_file = file(params.logo_file, checkIfExists: true)
css_file = file(params.css_file, checkIfExists: true)

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

include { TABULAR_TO_GSEA_CHIP } from '../modules/local/tabular_to_gsea_chip'

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
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main' 
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main' 
include { GSEA_GSEA                                         } from '../modules/nf-core/gsea/gsea/main' 
include { CUSTOM_TABULARTOGSEAGCT                           } from '../modules/nf-core/custom/tabulartogseagct/main'
include { CUSTOM_TABULARTOGSEACLS                           } from '../modules/nf-core/custom/tabulartogseacls/main' 
include { RMARKDOWNNOTEBOOK                                 } from '../modules/nf-core/rmarkdownnotebook/main' 

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
            if (!it.id){
                it.id = it.values().join('_')
            }
            it
        }

    // Firstly Filter the input matrix

    CUSTOM_MATRIXFILTER(
        VALIDATOR.out.fom.map{ tuple(it[0], it[3]) },
        VALIDATOR.out.fom.map{ tuple(it[0], it[1]) }
    )
    
    // Run the DESeq differential module, which doesn't take the feature
    // annotations 

    ch_samples_and_filtered_matrix = VALIDATOR.out.fom
        .map{
            tuple(it[0], it[1])                     // meta, samplesheet
        }
        .join(CUSTOM_MATRIXFILTER.out.filtered)     // -> meta, samplesheet, filtered matrix
        .map{ it.tail() } 

    DESEQ2_DIFFERENTIAL (
        ch_contrasts.combine(ch_samples_and_filtered_matrix),
        ch_control_features
    )

    // Run a gene set analysis where directed

    // Currently, we're letting GSEA work on the expression data. In future we
    // will allow use of GSEA preranked instead, which will work with the fold
    // changes/ p values from DESeq2
    
    if (params.gsea_run){    
    
        ch_gene_sets = Channel.from(gene_sets_file)

        // For GSEA, we need to convert normalised counts to a GCT format for
        // input, and process the sample sheet to generate class definitions
        // (CLS) for the variable used in each contrast        

        CUSTOM_TABULARTOGSEAGCT ( DESEQ2_DIFFERENTIAL.out.normalised_counts )

        ch_contrasts_and_samples = ch_contrasts.combine( VALIDATOR.out.fom.map { it[1] } )
        CUSTOM_TABULARTOGSEACLS(ch_contrasts_and_samples) 

        TABULAR_TO_GSEA_CHIP(
            VALIDATOR.out.fom.map{ it[2] },
            [params.features_id_col, params.features_name_col]    
        )

        ch_gsea_inputs = CUSTOM_TABULARTOGSEAGCT.out.gct
            .join(CUSTOM_TABULARTOGSEACLS.out.cls)
            .combine(ch_gene_sets)                        

        GSEA_GSEA( 
            ch_gsea_inputs,
            ch_gsea_inputs.map{ tuple(it[0].reference, it[0].target) }, // * 
            TABULAR_TO_GSEA_CHIP.out.chip.first()
        )
        
        // * Note: GSEA module currently uses a value channel for the mandatory
        // non-file arguments used to define contrasts, hence the indicated
        // usage of map to perform that transformation. An active subject of
        // debate

        ch_gsea_results = GSEA_GSEA.out.report_tsvs_ref
            .join(GSEA_GSEA.out.report_tsvs_target)

        // Record GSEA versions
        ch_versions = ch_versions
            .mix(TABULAR_TO_GSEA_CHIP.out.versions)
            .mix(GSEA_GSEA.out.versions)
    }

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

    ch_all_matrices = VALIDATOR.out.fom
        .combine(ch_processed_matrices)                         // Add processed marices to what we have in the FOM
        .map{
            tuple(it[0], it[1], it[2], [ it[3], it[4], it[5] ]) // Remove the experiment meta and group the matrices
        }
        .first()
 
    PLOT_EXPLORATORY(
        ch_contrast_variables
            .combine(ch_all_matrices.map{ it.tail() })
    )

    // Differential analysis using the results of DESeq2

    PLOT_DIFFERENTIAL(
        DESEQ2_DIFFERENTIAL.out.results, 
        ch_all_matrices
    )

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
    
    // Generate a list of files that will be used by the markdown report

    ch_report_file = Channel.from(report_file)
        .map{ tuple(exp_meta, it) }

    ch_logo_file = Channel.from(logo_file)
    ch_css_file = Channel.from(css_file)

    ch_report_input_files = ch_all_matrices
        .tail()
        .map{it.flatten()}
        .combine(ch_contrasts_file.map{it.tail()})
        .combine(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml)
        .combine(ch_logo_file)
        .combine(ch_css_file)
        .combine(DESEQ2_DIFFERENTIAL.out.results.map{it[1]}.toList())
 
    if (params.gsea_run){ 
        ch_report_input_files = ch_report_input_files
            .combine(ch_gsea_results
                .map{it.tail()}.flatMap().toList()
            )
    }

    // Make a params list - starting with the input matrices and the relevant
    // params to use in reporting

    ch_report_params = ch_report_input_files
        .map{[
            observations_file: it[0].name,
            features_file: it[1].name, 
            raw_matrix: it[2].name, 
            normalised_matrix: it[3].name, 
            variance_stabilised_matrix: it[4].name, 
            contrasts_file: it[5].name,
            versions_file: it[6].name,
            logo: it[7].name, 
            css: it[8].name
        ] + params.findAll{ k,v -> k.matches(~/^(study|observations|features|filtering|exploratory|differential|deseq2|gsea).*/) }}

    // TO DO: add further params - e.g. for custom logo etc, and for analysis
    // params 

    RMARKDOWNNOTEBOOK(
        ch_report_file,
        ch_report_params,
        ch_report_input_files
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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

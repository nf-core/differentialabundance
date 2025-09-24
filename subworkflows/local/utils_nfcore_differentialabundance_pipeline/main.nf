//
// Subworkflow with functionality specific to the nf-core/differentialabundance pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    // Check that params is available
    if (!params) {
        error("Pipeline parameters not initialized. This is a critical error.")
    }

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Get paramsets based on paramsheet or default parameters
    //
    paramsets = getParamsetConfigurations()
    ch_paramsets = Channel.fromList(paramsets)
        .map { paramset -> [
            id: paramset.study_name,
            paramset_name: paramset.paramset_name,
            params: paramset.findAll{ k,v -> k != 'paramset_name' }
        ]}

    //
    // Custom validate input parameters
    //
    validateInputParameters(paramsets)

    emit:
    paramsets = ch_paramsets
    versions  = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters(paramsets) {

    genomeExistsError()

    // check the params in each of the paramsets
    paramsets.each { row ->

        // check file existence for input parameters
        if (row.input) {
            file(row.input, checkIfExists: true)
        }

        // Validate study type spepecific parameters
        if (row.study_type == 'affy_array') {
            if (!row.affy_cel_files_archive) {
                error("CEL files archive not specified!")
            }
        } else if (row.study_type == 'maxquant') {
            if (row.functional_method) {
                error("Functional analysis is not yet possible with maxquant input data; please set --functional_method to null and rerun pipeline!")
            }
            if (!row.matrix) {
                error("Input matrix not specified!")
            }
        } else if (row.study_type == 'geo_soft_file') {
            if (!row.querygse || !row.features_metadata_cols) {
                error("Query GSE not specified or features metadata columns not specified")
            }
        } else if (row.study_type == "rnaseq") {
            if (!row.matrix) {
                error("Input matrix not specified!")
            }
        }

        // Validate functional analysis parameters
        if (row.functional_method) {
            if (row.functional_method == 'gsea' && !row.gene_sets_files) {
                error("GSEA activated but gene set file not specified!")
            } else if (row.functional_method == 'gprofiler2') {
                if (row.gene_sets_files) {
                    if (row.gene_sets_files.split(",").size() > 1) {
                        error("gprofiler2 can currently only work with a single gene set file")
                    }
                } else {
                    if (!row.gprofiler2_token && !row.gprofiler2_organism) {
                        error("To run gprofiler2, please provide a run token, GMT file or organism!")
                    }
                }
                // check if provided gprofiler2 background file exists
                if (row.gprofiler2_background_file && row.gprofiler2_background_file != 'auto') {
                    file(row.gprofiler2_background_file, checkIfExists: true)
                }
            }
        }

        // Validate contrasts parameters
        if (row.contrasts_yml && row.contrasts) {
            error("Both '--contrasts' and '--contrasts_yml' parameters are set. Please specify only one of these options to define contrasts.")
        }
        if (!(row.contrasts_yml || row.contrasts)) {
            error("Either '--contrasts' and '--contrasts_yml' must be set. Please specify one of these options to define contrasts.")
        }
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",

            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [

        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

// Get configurations based on whether use paramsheet or default params
def getParamsetConfigurations() {
    // Use paramsheet if paramset_name is provided, otherwise use default params
    def paramsets = (params.paramset_name) ? getParamsheetConfigurations() : getDefaultConfigurations()
    return paramsets.collect { paramset ->
        // some params are not useful through the pipeline run, remove them for cleaner meta
        def ignore = ['help', 'help_full', 'show_hidden', 'genomes']
        // non static params would interfere with cache
        // remove them from meta to avoid problems with resume
        def nonstaticparams = ['trace_report_suffix']
        paramset.findAll { k, v -> !(ignore + nonstaticparams).contains(k) } as Map
    }
}

// Get configurations from paramsheet
// To be able to retrieve a full set of parameters considering both paramsheet rows
// (highest priority) and pipeline params scope and validate them all together,
// we need to do the following:
// 1. get and validate paramsheet configurations
// 2. fill missing params with pipeline params
def getParamsheetConfigurations() {
    // get paramsheet path
    def paramsheet_path = file(params.paramsheet, checkIfExists: true)

    // Create temporary schema for validation - with only the fields in the paramsheet
    def schema_path = createParamsheetSchema(paramsheet_path)

    // Load paramsheet and validate each row against the transformed schema
    def raw_paramsheet = samplesheetToList(paramsheet_path, schema_path).collect { it -> it[0] }

    def paramsheet = raw_paramsheet
        // remove empty values
        .collect { row ->
            return row.findAll { key, value -> value != [] }
        }
        // Only keep row matching with paramset name
        .findAll { row ->
            return row.paramset_name == params.paramset_name
        }

    if (paramsheet.isEmpty()) {
        if (params.paramset_name) {
            error("No configuration found in paramsheet for paramset_name '${params.paramset_name}'")
        } else {
            error("No valid configurations found in paramsheet")
        }
    }

    // return paramset with paramsheet params (highest priority) and pipeline params
    return paramsheet
        .collect{ row ->
            // clean empty values with []
            def cleanparams = row.findAll { k, v -> v != [] } as Map
            // add missing params from pipeline params
            def fullparams = params + cleanparams
            return fullparams
        }
}

// Get default configurations from pipeline parameters
def getDefaultConfigurations() {
    // replace null by string 'contrasts' for paramset_name to avoid certain problems with null object
    return [params + [paramset_name: 'contrasts']]
}

// Create a temporary schema file for paramsheet validation
// This schema is derived from the pipeline's nextflow_schema.json by:
// 1. Reading the paramsheet headers to determine which fields to include
// 2. Extracting only those properties (also adding them as meta) and
// required fields from the pipeline schema
// 3. Creating a schema that allows for an array of objects
// @return The absolute path to the temporary schema file
def createParamsheetSchema(paramsheet_path) {
    // Load and parse pipeline schema
    def pipeline_schema = new File("${projectDir}/nextflow_schema.json").text
    def schema_json = new groovy.json.JsonSlurper().parseText(pipeline_schema)

    // Get headers from paramsheet
    def paramsheet_lines = paramsheet_path.newInputStream().withReader { it.readLines() }
    def headers = paramsheet_lines[0].split(',')
    // Ensure paramset_name is in headers
    if (!headers.contains('paramset_name')) {
        error("The paramsheet must contain a 'paramset_name' column")
    }

    // Extract all properties and required fields from schema
    def all_properties = extractPropertiesFromSchema(schema_json)
    def all_required = extractRequiredFromSchema(schema_json)

    // Restrict properties and required fields to those in the paramsheet
    def filtered_properties = all_properties.findAll { k, v -> headers.contains(k) }
    def filtered_required = all_required.findAll { k -> headers.contains(k) }

    // Create samplesheet schema with filtered properties and paramset_name as required
    def samplesheet_schema = [
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        'title': 'nf-core/differentialabundance - paramsheet schema',
        'description': 'Schema for validating the paramsheet configuration',
        'type': 'array',
        'items': [
            'type': 'object',
            'properties': filtered_properties,
            'required': filtered_required
        ]
    ]

    // Write temporary schema file
    def temp_schema = File.createTempFile("samplesheet_schema", ".json")
    def schema_string = new groovy.json.JsonBuilder(samplesheet_schema).toPrettyString()
    temp_schema.text = schema_string

    return temp_schema.absolutePath
}

// extract properties from schema
def extractPropertiesFromSchema(schema) {
    def all_properties = [:]
    schema.$defs.each { group_name, group_def ->
        if (group_def.properties) {
            group_def.properties.each { prop_name, prop_def ->
                // Add meta field to each property
                def prop_with_meta = prop_def + [meta: [prop_name]]
                all_properties[prop_name] = prop_with_meta
            }
        }
    }
    return all_properties
}

// extract required fields from schema
def extractRequiredFromSchema(schema) {
    def all_required = []
    schema.$defs.each { group_name, group_def ->
        if (group_def.required) {
            group_def.required.each { req_prop ->
                all_required << req_prop
            }
        }
    }
    if (!all_required.contains('paramset_name')) {
        all_required << 'paramset_name'
    }
    return all_required
}

// prepare the input for the module by keeping only the relevant params
// in the meta and group channel by simplified meta and unique inputs.
// By keeping only the relevant params in meta when calling the modules,
// we could ensure proper functionality of pipeline `resume`.
// @param channel: the input channel
// @param category: the category in which the module belong to
// @return a channel with [simplified meta, files ...]
// the simplified meta would have the following structure:
// [paramset_names: [list of paramset names], params: [relevant params map], ...other k,v]
def prepareModuleInput(channel, category) {
    return channel
        .map {
            // clean the params by keeping only the relevant params
            def simplifiedparams = getRelevantParams(it[0].params, category)
            // replace meta.params by simplified params
            def simplifiedmeta = it[0] + [params: simplifiedparams]

            // use simplified meta as key
            [simplifiedmeta, it[0].paramset_name, it[1..-1]]  // [ meta, paramset_name, [files] ]
        }
        .groupTuple()
        .map { simplifiedmeta, paramset_names, file_lists ->
            // replace paramset_name by paramset_names
            def meta = [paramset_names: paramset_names] + simplifiedmeta.findAll{ k,v -> k != 'paramset_name'}
            // the list files are the same for the same simplified meta,
            // thus the list of files generated from grouping are just a repetition of the same files
        [meta] + file_lists[0]
        }
}

// prepare the output for the module by adding back the full paramsets
// to the meta. This is done by matching with paramset_name.
// When provided, it will also remove the unnecesary keys from meta.
// @param channel: the output channel from module
// @param paramsets: the channel containing all the parameter sets
// @param meta_keys_to_remove: a list of meta keys to remove. null if not
// provided.
// @param use_meta_key: use the base meta structure as key. So
// [id: study_name, paramset_name: paramset_name, params: [:]]
def prepareModuleOutput(channel, paramsets, List meta_keys_to_remove = null, Boolean use_meta_key = false) {
    // Prepare the channel to have paramset_name as key
    channel_by_name = channel
        // First parse the channel to have paramset_names as key
        .map { [it[0].paramset_names] + it }
        // Transpose the list of paramset_names
        .transpose(by:0)
        .map {
            def paramset_name = it[0]
            // Replace paramset_names by paramset_name in the meta
            def meta = [paramset_name: paramset_name] + it[1].findAll{ k,v -> k != 'paramset_names' }   // [paramset_name, params, ...]
            // Put paramset name as key
            [paramset_name, meta] + it[2..-1]  // [paramset_name, meta, files ...]
        }
    // Add back the full paramset by matching with paramsets channel with paramset_name as key
    channel_with_updated_paramsets = paramsets
        .map { [it.paramset_name, it] }
        // We use combine instead of join, because for the same paramset_name,
        // different meta and files may have been generated because of the different contrasts
        .combine(channel_by_name, by:0)
        .map {
            def meta_paramset = it[1]
            def meta_out = it[2]
            // Remove unnecessary keys from meta, when asked
            def meta_cleaned = (meta_keys_to_remove) ? meta_out.findAll{ k,v -> !meta_keys_to_remove.contains(k) } : meta_out
            // Replace output meta simplified params by full params from paramset
            def meta = meta_cleaned + [params: meta_paramset.params]

            if (use_meta_key) {
                // Define a key using the basic meta structure: only containing id, paramset_name and params, when asked.
                // Note that all the channels in the pipeline have study_name as id, except those containing contrast info.
                // Hence, we need to use the study_name as id in the key.
                def key = [id: meta.params.study_name, paramset_name: meta.paramset_name, params: meta.params]
                [key, meta] + it[3..-1] // [key, meta with full paramset, files ...]
            } else {
                [meta] + it[3..-1]      // [meta with full paramset, files ...]
            }
        }
}

// get the relevant params given the category
// This function will parse the schema to only retain the parameters that are relevant.
// The relevant parameters include: 1) those parameters from the base category that
// might be shared by many modules in the pipeline, 2) those parameters needed for
// preceding categories and 3) those parameters that are only needed for a given
// category. The category is defined in the parameter group name in the schema.
// @param meta: the meta map
// @param category: the category name
def getRelevantParams(paramset, category) {
    // Define schema URL - in practice this would be loaded from file
    def schema = new groovy.json.JsonSlurper().parseText(new File("${projectDir}/nextflow_schema.json").text)

    // the relevant groups are the one defined by the category
    // and all the preceding ones
    def relevant_categories = [
        'preprocessing': ['base', 'preprocessing'],
        'exploratory': ['base', 'preprocessing', 'exploratory'],
        'differential': ['base', 'preprocessing', 'differential'],
        'functional': ['base', 'preprocessing', 'differential', 'functional']
    ]
    if (!relevant_categories.containsKey(category)) {
        error("Category '${category}' not found in schema.")
    }

    // Get all relevant categories for this level
    def categories_to_check = relevant_categories[category]

    // Initialize results map
    def relevantParams = [:]

    // For each relevant category, find and collect parameters
    categories_to_check.each { cat ->
        def matchingGroups = schema['$defs'].findAll { groupName, group ->
            groupName.startsWith("${cat}_")
        }

        matchingGroups.each { groupName, group ->
            if (group.properties) {
                group.properties.each { paramName, paramProps ->
                    if (paramset.containsKey(paramName)) {
                        relevantParams[paramName] = paramset[paramName]
                    }
                }
            }
        }
    }

    return relevantParams
}

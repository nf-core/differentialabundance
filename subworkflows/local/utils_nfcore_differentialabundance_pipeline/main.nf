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
include { validate                  } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
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
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

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

    def before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/differentialabundance ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    def after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x
* Software dependencies
    https://github.com/nf-core/differentialabundance/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll('\\u001b\\[[0-9;]*m', '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        false,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
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
    def configurations = params.paramsheet
        ? getParamsheetConfigurations()
        : getDefaultConfigurations()
    paramsets = validateConfigurations(configurations)
        .collect { paramset -> addDifferentialRuntimeParams(paramset) }
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

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Per-method DISPLAY columns stamped into meta.params at paramset-parse
// time, used downstream by report.qmd, shinyngs and the output DSL.
// FILTER columns are owned separately by the abundance_differential_filter
// subworkflow's getDifferentialMethodParams (deliberately different for
// propd: `significant` is the filter column, `rcDdis` is the display score).
def getDifferentialMethodRuntimeParams(differential_method) {
    def runtime_params = [
        'deseq2': [
            differential_fc_column         : 'log2FoldChange',
            differential_pval_column       : 'pvalue',
            differential_qval_column       : 'padj',
            differential_foldchanges_logged: true
        ],
        'limma' : [
            differential_fc_column         : 'logFC',
            differential_pval_column       : 'P.Value',
            differential_qval_column       : 'adj.P.Val',
            differential_foldchanges_logged: true
        ],
        'propd' : [
            differential_fc_column         : 'LFC',
            differential_pval_column       : 'rcDdis',
            differential_qval_column       : 'rcDdis',
            differential_foldchanges_logged: true
        ],
        'dream' : [
            differential_fc_column         : 'logFC',
            differential_pval_column       : 'P.Value',
            differential_qval_column       : 'adj.P.Val',
            differential_foldchanges_logged: true
        ]
    ][differential_method]

    runtime_params ?: [:]
}

def addDifferentialRuntimeParams(paramset) {
    def runtime_params = getDifferentialMethodRuntimeParams(paramset.differential_method)
    def missing_runtime_params = runtime_params.findAll { key, _value ->
        !paramset.containsKey(key)
    }
    paramset + missing_runtime_params
}

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
                error("CEL files archive not specified for paramset={${row.paramset_name}}!")
            }
        } else if (row.study_type == 'maxquant') {
            if (row.functional_method && row.functional_method != 'none') {
                error("Functional analysis is not yet possible with maxquant input data; please set --functional_method to 'none' and rerun pipeline!")
            }
            if (!row.matrix) {
                error("Input matrix not specified for paramset={${row.paramset_name}}!")
            }
        } else if (row.study_type == 'geo_soft_file') {
            if (!row.querygse || !row.features_metadata_cols) {
                error("Query GSE not specified or features metadata columns not specified for paramset={${row.paramset_name}}")
            }
        } else if (row.study_type in ["rnaseq", "generic_matrix"]) {
            if (!row.matrix) {
                error("Input matrix not specified for paramset={${row.paramset_name}}!")
            }
        }

        // Validate functional analysis parameters
        if (row.functional_method && row.functional_method != 'none') {
            if (row.functional_method == 'gsea' && !row.gene_sets_files) {
                error("GSEA activated but gene set file not specified for paramset={${row.paramset_name}}!")
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
            } else if (row.functional_method == 'decoupler') {
                if (!row.decoupler_network) {
                        error("To run decoupler, please provide a network file!")
                    }
            }
        }

        // Validate contrasts parameters
        if (!row.contrasts) {
            error("'--contrasts' must be set. Please provide a contrasts file in CSV, TSV, YML, or YAML format.")
        }

        // Validate control-based size factor parameters
        if (row.sizefactors_from_controls && row.differential_method != 'deseq2') {
            error("'--sizefactors_from_controls' is only supported with '--differential_method deseq2'. Found differential_method='${row.differential_method}' in paramset='${row.paramset_name}'.")
        }

        if (row.control_features && !row.sizefactors_from_controls && row.differential_method != 'deseq2') {
            log.warn("'--control_features' with '${row.differential_method}' will only strip control features from matrices. Normalization from controls is only supported with '--differential_method deseq2'.")
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

// Validate configurations against the schema.
def validateConfigurations(configurations) {
    return configurations.collect { paramset ->
        // Some params are not useful through the pipeline run.
        // Remove them for cleaner meta
        def ignore = ['help', 'help_full', 'show_hidden', 'genomes']
        // Non static params would interfere with cache.
        // Remove them from meta to avoid problems with resume
        def nonstaticparams = ['trace_report_suffix']
        def cleanparamset = paramset.findAll { k, v -> !(ignore + nonstaticparams).contains(k) } as Map

        // Remove null to skip validation on them
        // This is needed because validate() will fail otherwise
        def notnullparams = cleanparamset.findAll { k, v -> v != null } as Map

        try {
            // Validate against schema
            validate(notnullparams, "${projectDir}/nextflow_schema.json")
        } catch (e) {
            // Surface the paramset name; nf-schema will then produce a detailed error.
            log.error "Validation failed for paramsheet row: ${paramset.paramset_name}"
            validate(notnullparams, "${projectDir}/nextflow_schema.json")
        }

        return cleanparamset
    }
}

// Get configurations from paramsheet
def getParamsheetConfigurations() {
    // Get paramsheet path
    def paramsheet_path = file(params.paramsheet, checkIfExists: true)

    // Load configs from paramsheet as list of maps
    def raw_paramsheet = loadYamlConfigs(paramsheet_path)

    def paramsheet = raw_paramsheet
        .findAll { row ->
            if (!params.paramset_name || params.paramset_name == 'all') {
                return true
            } else {
                // Only keep row matching with paramset name(s)
                return row.paramset_name in params.paramset_name.tokenize(',')
            }
        }

    if (paramsheet.isEmpty()) {
        if (params.paramset_name) {
            error("No configuration found in paramsheet for paramset_name '${params.paramset_name}'")
        } else {
            error("No valid configurations found in paramsheet")
        }
    }

    return paramsheet
        .collect{ row ->
            // Note that the paramsheet may not contain all the parameters
            // defined in the pipeline, so we need to merge them
            def fullparamset = params + row
            return fullparamset
        }
}

// Get default configurations from pipeline parameters (profile mode)
def getDefaultConfigurations() {
    // Use paramset_name from profile if set, otherwise fall back to 'contrasts'
    def pname = params.paramset_name ?: 'contrasts'
    return [params + [paramset_name: pname]]
}

// Load configurations from yaml file
def loadYamlConfigs(yaml_path) {
    // Load yaml content
    def configs = loadYaml(yaml_path)

    // Resolve includes for each config
    configs = configs.collect { config ->
        config = resolveIncludes(config)
        config.remove('include')
        return config
    }

    return configs
}
def loadYaml(yaml_path) {
    // Load yaml content
    def yaml_content = yaml_path.text

    // Substitute ${projectDir} with actual value
    // alternative ways? This can be fragile
    yaml_content = yaml_content.replaceAll('\\$\\{projectDir\\}', projectDir.toString())
    yaml_content = yaml_content.replaceAll('\\$projectDir', projectDir.toString())

    // Parse yaml content
    def yaml_parser = new org.yaml.snakeyaml.Yaml()
    def loaded = yaml_parser.load(yaml_content)
    def configs = (loaded instanceof List) ? loaded : [loaded]

    return configs
}
// Helper function to resolve includes recursively
def resolveIncludes(config) {
    def yaml_parser = new org.yaml.snakeyaml.Yaml()

    if (config.containsKey('include')) {
        def includePaths = config.include.split(',') // Split by comma to handle multiple includes
        includePaths.each { includePath ->
            def includeParts = includePath.split('/')
            def includeFile = includeParts[0]
            def paramsetName = includeParts[1]

            // Load the included YAML file
            def includeFilePath = file("${projectDir}/conf/${includeFile}.yaml")
            if (!includeFilePath.exists()) {
                error("Included file '${includeFilePath}' not found.")
            }
            def includeConfigs = loadYaml(includeFilePath)

            // Find the paramset with the given name
            def includedConfig = includeConfigs.find { it.paramset_name == paramsetName }
            if (!includedConfig) {
                error("Paramset '${paramsetName}' not found in included file '${includeFilePath}'.")
            }

            // Recursively resolve includes in the included config
            includedConfig = resolveIncludes(includedConfig)

            // Merge configs
            includedConfig.putAll(config)
            config = includedConfig
        }
    }
    return config
}

// prepare the input for the module by keeping only the relevant params
// in the meta and group channel by simplified meta and unique inputs.
// By keeping only the relevant params in meta when calling the modules,
// we could ensure proper functionality of pipeline `resume`.
// @param channel: the input channel
// @param category: the category in which the module belong to
// @return a channel with [simplified meta, files ...]
// the simplified meta would have the following structure:
// [params: [relevant params map], ...other k,v]
def prepareModuleInput(channel, category) {
    return channel
        .map {
            // clean the params by keeping only the relevant params
            def simplifiedparams = getRelevantParams(it[0].params, category)
            // replace meta.params by simplified params
            def simplifiedmeta = it[0] + [params: simplifiedparams]

            // remove paramset_name from meta, and use as key
            // This avoids the paramset_name id being considered as a parameter that can interfere with caching and resume
            def key = simplifiedmeta.findAll { k, v -> k != 'paramset_name' }
            [key, it[1..-1]]  // [ key, [files] ]
        }
        .groupTuple()
        .map { simplifiedmeta, file_lists ->
            // the list files are the same for the same simplified meta,
            // thus the list of files generated from grouping are just a repetition of the same files
            [simplifiedmeta] + file_lists[0]
        }
}

// prepare the output for the module by adding back the full paramsets
// to the meta. This is done by making sure all the k,v from the
// channel meta match with those from the full paramset; if so we add
// back the rest of params from the full paramset.
// @param channel: the output channel from module
// @param paramsets: the channel containing all the parameter sets
// @param meta_keys_to_remove: (optional) a list of meta keys to remove.
// @param use_meta_key: (optional) use the base meta structure as key. So
// [id: study_name, paramset_name: paramset_name, params: [:]] without other
// k,v from the meta.
def prepareModuleOutput(channel, paramsets, List meta_keys_to_remove = null, Boolean use_meta_key = false) {
    return paramsets
        .combine(channel)  // [ meta_paramset, meta_out, files ... ]
        .filter { it ->
            def meta_paramset = it[0]
            def meta_out = it[1]
            // check if all the k,v in meta_out.params match with those in meta_paramset.params
            meta_out.params.every { key, value -> meta_paramset.params[key] == value }
        }
        .map { it ->
            def meta_paramset = it[0]
            def meta_out = it[1]

            // Remove unnecessary keys from meta, when asked
            def meta_cleaned = meta_keys_to_remove ? meta_out.findAll { k, v -> !meta_keys_to_remove.contains(k) } : meta_out

            // Replace output meta simplified params by full params from paramset
            def meta = meta_cleaned + [paramset_name: meta_paramset.paramset_name, params: meta_paramset.params]

            if (use_meta_key) {
                // Define a key using the basic meta structure: only containing id, paramset_name and params, when asked.
                // Note that all the channels in the pipeline have study_name as id, except those containing contrast info.
                // Hence, we need to use the study_name as id in the key.
                def key = [id: meta.params.study_name, paramset_name: meta.paramset_name, params: meta.params]
                [key, meta] + it[2..-1] // [key, meta with full paramset, files ...]
            } else {
                [meta] + it[2..-1]      // [meta with full paramset, files ...]
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

    [
        'differential_fc_column',
        'differential_pval_column',
        'differential_qval_column',
        'differential_foldchanges_logged'
    ].each { paramName ->
        if (paramset.containsKey(paramName)) {
            relevantParams[paramName] = paramset[paramName]
        }
    }

    return relevantParams
}

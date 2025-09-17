/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { prepareModuleInput     } from '../subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { prepareModuleOutput    } from '../subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { UNTAR                                             } from '../modules/nf-core/untar/main.nf'
include { SHINYNGS_APP                                      } from '../modules/nf-core/shinyngs/app/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../modules/nf-core/shinyngs/staticexploratory/main'
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../modules/nf-core/shinyngs/staticdifferential/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main'
include { RMARKDOWNNOTEBOOK                                 } from '../modules/nf-core/rmarkdownnotebook/main'
include { AFFY_JUSTRMA as AFFY_JUSTRMA_RAW                  } from '../modules/nf-core/affy/justrma/main'
include { AFFY_JUSTRMA as AFFY_JUSTRMA_NORM                 } from '../modules/nf-core/affy/justrma/main'
include { PROTEUS_READPROTEINGROUPS as PROTEUS              } from '../modules/nf-core/proteus/readproteingroups/main'
include { GEOQUERY_GETGEO                                   } from '../modules/nf-core/geoquery/getgeo/main'
include { ZIP as MAKE_REPORT_BUNDLE                         } from '../modules/nf-core/zip/main'
include { IMMUNEDECONV                                      } from '../modules/nf-core/immunedeconv/main'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// SUBWORKFLOW: Installed directly from nf-core/modules
//
include { ABUNDANCE_DIFFERENTIAL_FILTER                     } from '../subworkflows/nf-core/abundance_differential_filter/main'
include { DIFFERENTIAL_FUNCTIONAL_ENRICHMENT                } from '../subworkflows/nf-core/differential_functional_enrichment/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    take:
    ch_paramsets

    main:

    ch_versions = Channel.empty()

    // ========================================================================
    // Handle input
    // ========================================================================

    // Get the sample sheets
    ch_samplesheet = ch_paramsets
        .map { meta ->
            [ meta, file(meta.params.input, checkIfExists: true) ]
        }

    // Create input channels based on study type
    ch_input = ch_samplesheet
        .branch { meta, input ->
            affy_array: meta.params.study_type == 'affy_array'
            maxquant: meta.params.study_type == 'maxquant'
            geo_soft_file: meta.params.study_type == 'geo_soft_file'
            rnaseq: meta.params.study_type == 'rnaseq'
        }

    // Handle Affy array inputs
    ch_celfiles = ch_input.affy_array
        .map { meta, input ->
            [ meta, file(meta.params.affy_cel_files_archive, checkIfExists: true) ]
        }

    // Handle Maxquant inputs
    ch_maxquant_inputs = ch_input.maxquant
        .map { meta, input ->
            [ meta, input, file(meta.params.matrix, checkIfExists: true) ]
        }

    // Handle GEO soft file inputs
    ch_querygse = ch_input.geo_soft_file
        .map { meta, input ->
            [ meta, meta.params.querygse ]
        }

    // Create optional parameter channels based on ch_paramsets
    ch_transcript_lengths = ch_paramsets
        .map { meta ->
            [ meta, meta.params.transcript_length_matrix ? file(meta.params.transcript_length_matrix, checkIfExists: true) : [] ]
        }

    ch_control_features = ch_paramsets
        .map { meta ->
            [ meta, meta.params.control_features ? file(meta.params.control_features, checkIfExists: true) : [] ]
        }

    ch_gene_sets = ch_paramsets
        .map { meta ->
            if (meta.params.functional_method == 'decoupler' && meta.params.decoupler_network) {
                [ meta, [file(meta.params.decoupler_network, checkIfExists: true)] ]
            } else {
                [ meta, meta.params.gene_sets_files ? meta.params.gene_sets_files.split(",").collect { file(it, checkIfExists: true) } : [] ]
            }
        }

    // ========================================================================
    // Handle contrasts
    // ========================================================================

    // Create contrasts channels
    ch_contrasts_file = ch_paramsets
        .map { meta ->
            [ meta, file(meta.params.contrasts_yml ?: meta.params.contrasts, checkIfExists: true) ]
        }

    ch_contrasts_file_with_extension = ch_contrasts_file
        .map { meta, file ->
            [ meta, file, file.extension ]
        }

    ch_contrast_variables_input = ch_contrasts_file_with_extension
        .branch{ meta, file, extension ->
            yml: extension == 'yml' || extension == 'yaml'
            csv: extension == 'csv'
            tsv: extension == 'tsv'
        }

    ch_contrasts_variables_from_yml = ch_contrast_variables_input.yml
        .flatMap { meta, yaml_file, ext ->
            def yaml_data = new groovy.yaml.YamlSlurper().parse(yaml_file)
            yaml_data.contrasts.collect { contrast ->
                if (contrast.containsKey('formula')) {
                return null
                }
                else if (contrast.containsKey('comparison')) { //  Necessary line for Maxquant to work. Check if it can be simplified to use contrast.id
                    tuple(meta, contrast.comparison[0])
                }
            }
        }
            .filter { it != null }
        .unique() // Uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)

    ch_contrasts_variables_from_other = ch_contrast_variables_input.csv.splitCsv(header:true)
        .mix(ch_contrast_variables_input.tsv.splitCsv(header:true, sep:'\t'))
        .map { meta, row, ext ->
            [ meta, row.variable ]
        }
        .unique()

    ch_contrast_variables = ch_contrasts_variables_from_yml
        .mix(ch_contrasts_variables_from_other)

    // ========================================================================
    // Data type specific preprocessing
    // ========================================================================

    //
    // 1. Deal with affy array data
    // If we have affy array data in the form of CEL files we'll be deriving
    // matrix and annotation from them
    //

    // Uncompress the CEL files archive
    UNTAR ( prepareModuleInput(ch_celfiles, 'preprocessing') )
    ch_untar_out = prepareModuleOutput(UNTAR.out.untar, ch_paramsets)

    // Run affy

    ch_affy_input = prepareModuleInput(ch_input.affy_array.join(ch_untar_out), 'preprocessing')

    AFFY_JUSTRMA_RAW (
        ch_affy_input,
        [[],[]]
    )
    AFFY_JUSTRMA_NORM (
        ch_affy_input,
        [[],[]]
    )

    ch_affy_raw = prepareModuleOutput(AFFY_JUSTRMA_RAW.out.expression, ch_paramsets)
    ch_affy_norm = prepareModuleOutput(AFFY_JUSTRMA_NORM.out.expression, ch_paramsets)
    ch_affy_platform_features = prepareModuleOutput(AFFY_JUSTRMA_RAW.out.annotation, ch_paramsets)

    ch_versions = ch_versions
        .mix(AFFY_JUSTRMA_RAW.out.versions)

    //
    // 2. Deal with maxquant data
    // We'll be running Proteus once per unique contrast variable to generate plots
    //

    // Add contrast variable
    ch_proteus_input = ch_maxquant_inputs
        .combine(ch_contrast_variables, by:0)
        .map { meta, input, matrix, contrast ->
            def meta_new = meta + [contrast: contrast]
            [meta_new, input, matrix]
        }

    // Run proteus to import protein abundances
    PROTEUS( prepareModuleInput(ch_proteus_input, 'preprocessing') )

    // Note that the tables are the same across contrasts, only one table will be necessary
    // that is why here we take the first one and remove the contrast variable from meta
    ch_proteus_raw = prepareModuleOutput(PROTEUS.out.raw_tab, ch_paramsets, meta_keys_to_remove=['contrast'])
        .first()
    ch_proteus_norm = prepareModuleOutput(PROTEUS.out.norm_tab, ch_paramsets, meta_keys_to_remove=['contrast'])
        .first()

    ch_versions = ch_versions.mix(PROTEUS.out.versions)

    //
    // 3. Deal with GEO soft file data
    // Run GEO query to get the annotation
    //

    GEOQUERY_GETGEO( prepareModuleInput(ch_querygse, 'preprocessing') )

    ch_soft_norm = prepareModuleOutput(GEOQUERY_GETGEO.out.expression, ch_paramsets)
    ch_soft_features = prepareModuleOutput(GEOQUERY_GETGEO.out.annotation, ch_paramsets)

    ch_versions = ch_versions
        .mix(GEOQUERY_GETGEO.out.versions)

    // ========================================================================
    // Parse channels for input matrices and features
    // ========================================================================

    //
    // Define raw and normalised input expression/abundance data
    //

    // Raw inputs

    ch_in_raw = ch_input.rnaseq
        .map{meta, input -> [meta, file(meta.params.matrix, checkIfExists: true)]}
        .mix(ch_affy_raw)
        .mix(ch_proteus_raw)

    // Normalised inputs

    ch_in_norm = ch_affy_norm
        .mix(ch_proteus_norm)
        .mix(ch_soft_norm)

    //
    // Fetch or derive a feature annotation table
    //

    // Branch ch_paramsets based on feature source
    ch_feature_sources = ch_paramsets
        .branch {
            user_features: it.params.features
            affy_features: it.params.study_type == 'affy_array'
            geo_features: it.params.study_type == 'geo_soft_file'
            gtf_features: it.params.gtf
            matrix_features: true
        }

    // Handle user-provided feature annotations
    ch_user_features = ch_feature_sources.user_features
        .map { meta -> [ meta, file(meta.params.features, checkIfExists: true) ] }

    // Handle Affy array platform features
    ch_affy_features = ch_feature_sources.affy_features
        .join(ch_affy_platform_features)

    // Handle GEO soft file features
    ch_geo_features = ch_feature_sources.geo_features
        .join(ch_soft_features)

    // Handle GTF-based feature annotations
    ch_gtf_files = ch_feature_sources.gtf_features
        .map { meta -> [ meta, file(meta.params.gtf, checkIfExists: true) ] }

    // Process GTF files if needed
    ch_gtf_for_processing = ch_gtf_files
        .branch {
            compressed: it[1].name.endsWith('.gz')
            uncompressed: !it[1].name.endsWith('.gz')
        }

    // Decompress GTF files if needed
    GUNZIP_GTF( prepareModuleInput(ch_gtf_for_processing.compressed, 'preprocessing') )
    ch_gunzip_out = prepareModuleOutput(GUNZIP_GTF.out.gunzip, ch_paramsets)
    ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

    // Combine compressed and uncompressed GTF files
    ch_gtf_processed = ch_gtf_for_processing.uncompressed
        .mix(ch_gunzip_out)

    // Convert GTF to feature annotation table
    GTF_TO_TABLE(
        prepareModuleInput(ch_gtf_processed, 'preprocessing'),
        [tuple('id':""), []]
    )
    ch_gtf_features = prepareModuleOutput(GTF_TO_TABLE.out.feature_annotation, ch_paramsets)
    ch_versions = ch_versions.mix(GTF_TO_TABLE.out.versions)

    // Extract features from matrix
    // Note that in the case of maxquant we use the normalised matrix
    // Whereas for other study types we use the raw matrix
    ch_pre_matrix_features = ch_feature_sources.matrix_features.branch{ meta ->
        maxquant: meta.params.study_type == 'maxquant'
        other: true
    }
    ch_matrix_features = ch_pre_matrix_features.maxquant.join(ch_in_norm)
        .mix(ch_pre_matrix_features.other.join(ch_in_raw))
        .map{ meta, matrix ->
            def matrix_as_anno_filename = "${workflow.workDir}/${matrix.getBaseName()}_as_anno.${matrix.getExtension()}"
            def matrix_copy = file(matrix_as_anno_filename)
            matrix_copy.exists() && matrix.getText().md5().equals(matrix_copy.getText().md5()) ?: matrix.copyTo(matrix_as_anno_filename)
            [ meta, file(matrix_as_anno_filename) ]
        }

    // Combine all feature annotation channels
    ch_features = ch_user_features
        .mix(ch_affy_features)
        .mix(ch_geo_features)
        .mix(ch_gtf_features)
        .mix(ch_matrix_features)

    // ========================================================================
    // Validate input
    // ========================================================================

    // Check compatibility of FOM elements and contrasts

    ch_matrix_sources = ch_paramsets.branch{ meta ->
        affy_or_maxquant: meta.params.study_type == 'affy_array' || meta.params.study_type == 'maxquant'
        geo_soft_file: meta.params.study_type == 'geo_soft_file'
        other: true
    }

    // Create a channel with the matrices for validation. For Affy and MaxQuant
    // we have both raw and norm matrices, for GEO soft file we only have norm,
    // for other we only have raw

    ch_matrices_for_validation = ch_matrix_sources.affy_or_maxquant
        .join(ch_in_raw)
        .join(ch_in_norm)
        .map{tuple(it[0], [it[1], it[2]])}
        .mix(ch_matrix_sources.geo_soft_file.join(ch_in_norm))
        .mix(ch_matrix_sources.other.join(ch_in_raw))

    // Validate the components against one another. We try not to assume too
    // much about chanel order, so we deploy the multiMap workaround.

    validator_input = ch_samplesheet
        .join(ch_matrices_for_validation)
        .join(ch_features)
        .join(ch_contrasts_file)

    validator_input = prepareModuleInput(validator_input, 'preprocessing')
        .multiMap{meta, samplesheet, matrices, features_file, contrasts_file ->
            samplesheet_matrices: [meta, samplesheet, matrices]
            features_file: [meta, features_file]
            contrasts_file: [meta, contrasts_file]
        }

    VALIDATOR(
        validator_input.samplesheet_matrices,
        validator_input.features_file,
        validator_input.contrasts_file
    )

    ch_validated_assays = prepareModuleOutput(VALIDATOR.out.assays, ch_paramsets)
    ch_validated_contrast = prepareModuleOutput(VALIDATOR.out.contrasts, ch_paramsets)
    ch_validated_samplemeta = prepareModuleOutput(VALIDATOR.out.sample_meta, ch_paramsets)
    ch_validated_featuremeta = prepareModuleOutput(VALIDATOR.out.feature_meta, ch_paramsets)

    // For Affy, we've validated multiple input matrices for raw and norm,
    // we'll separate them out again here

    ch_multi_validated_assays = ch_validated_assays
        .filter{meta, assays -> meta.params.study_type == 'affy_array' || meta.params.study_type == 'maxquant'}
        .transpose()
        .branch { meta, assay ->
            raw: assay.name.contains('raw')
            normalised: assay.name =~ /normali[sz]ed/
        }

    // Get raw matrices from the validation

    ch_raw = ch_multi_validated_assays.raw
        .mix(ch_validated_assays.filter{meta, assay -> meta.params.study_type == 'rnaseq'})

    // For RNASeq and GEO soft file we've validated a single matrix, raw in the
    // case of RNASeq and norm in the case of GEO soft file, and these are the
    // matrices we'll use for differential abundance testing. For affy or maxquant
    // we'll use the normalised matrices.

    ch_matrix_for_differential = ch_multi_validated_assays.normalised
        .mix(ch_validated_assays.filter{meta, assay -> meta.params.study_type == 'rnaseq' || meta.params.study_type == 'geo_soft_file'})

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = ch_validated_contrast
        .splitCsv ( header:true, sep:'\t' )
        .map{meta, contrast ->
            contrast.blocking = contrast.blocking.replaceAll('^NA$', '')
            if (!contrast.id){
                contrast.id = contrast.values().join('_')
            }
            contrast.formula = contrast.formula?.trim() ? contrast.formula.trim() : null
            contrast.make_contrasts_str = contrast.make_contrasts_str?.trim() ? contrast.make_contrasts_str.trim() : null

            return [meta, contrast, contrast.variable, contrast.reference, contrast.target, contrast.formula, contrast.make_contrasts_str]
        }
        .groupTuple() // [meta, [contrast], [variable], [reference], [target], [formula], [comparison]]

    // =======================================================================
    // Single cell deconvolution
    // =======================================================================

    ch_immunedeconv_input = ch_in_raw
        .filter{meta, raw -> meta.params.immunedeconv_run}

    ch_immunedeconv_input = prepareModuleInput(ch_immunedeconv_input, 'preprocessing')
        .multiMap{meta, raw ->
            input: [meta, raw, meta.params.immunedeconv_method, meta.params.immunedeconv_function]
            name_col: meta.params.features_name_col
        }

    IMMUNEDECONV(
        ch_immunedeconv_input.input,
        ch_immunedeconv_input.name_col
    )

    ch_versions = ch_versions
        .mix(IMMUNEDECONV.out.versions)

    // ========================================================================
    // Filter matrix
    // ========================================================================

    ch_matrixfilter_input = ch_matrix_for_differential
        .join(ch_validated_samplemeta)

    ch_matrixfilter_input = prepareModuleInput(ch_matrixfilter_input, 'preprocessing')
        .multiMap{meta, matrix, samplesheet ->
            matrix: [meta, matrix]
            samplesheet: [meta, samplesheet]
        }

    CUSTOM_MATRIXFILTER(
        ch_matrixfilter_input.matrix,
        ch_matrixfilter_input.samplesheet
    )

    ch_filtered_matrix = prepareModuleOutput(CUSTOM_MATRIXFILTER.out.filtered, ch_paramsets)

    // ========================================================================
    // Differential analysis
    // ========================================================================

    ch_differential_input = ch_filtered_matrix
        .join(ch_validated_samplemeta)
        .join(ch_transcript_lengths)
        .join(ch_control_features)
        .join(ch_contrasts)

    // Use a multiMap to generate synched channels for differential analysis
    ch_differential_input = prepareModuleInput(ch_differential_input, 'differential')
        .multiMap{ meta, matrix, samplesheet, transcript_lengths, control_features, contrast, variable, reference, target, formula, comparison ->
            input: [meta, matrix, meta.params.differential_method, meta.params.differential_min_fold_change, meta.params.differential_max_qval]
            samplesheet: [meta, samplesheet]
            transcript_lengths: [meta, transcript_lengths]
            control_features: [meta, control_features]
            contrast: [meta, contrast, variable, reference, target, formula, comparison]
        }

    // Run differential analysis

    ABUNDANCE_DIFFERENTIAL_FILTER(
        ch_differential_input.input,
        ch_differential_input.samplesheet,
        ch_differential_input.transcript_lengths,
        ch_differential_input.control_features,
        ch_differential_input.contrast
    )

    // Collect differential results
    // Note that 'differential_method' is additionally added to the meta in the differential subworkflow.
    // Remove it to keep a consistent meta structure (needed for join/combine).

    // Also note that these channels, the meta contain other info apart from the base paramset meta.
    // Hence we create a key using the simple paramset meta with only meta.id, meta.paramset_name and meta.params,
    // by setting 'use_meta_key' to true. This will facilitate later on to join/combine channels.
    ch_differential_results = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise, ch_paramsets, meta_keys_to_remove=['differential_method'], use_meta_key=true) // key, meta, results
    ch_differential_results_filtered = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise_filtered, ch_paramsets, meta_keys_to_remove=['differential_method'], use_meta_key=true) // key, meta, results_filtered
    ch_differential_model = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.model, ch_paramsets, meta_keys_to_remove=['differential_method'], use_meta_key=true) // key, meta, model

    // Whereas these channels, the meta do not contain contrast info, as they come from the NORM modules instead of DIFFERENTIAL modules.
    // We do not need to define an extra key based on meta for later join/combine.
    ch_differential_norm = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.normalised_matrix, ch_paramsets, meta_keys_to_remove=['differential_method']) // meta, norm file
    ch_differential_varstab = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.variance_stabilised_matrix, ch_paramsets, meta_keys_to_remove=['differential_method']) // meta, varstab file

    ch_versions = ch_versions
        .mix(ABUNDANCE_DIFFERENTIAL_FILTER.out.versions)

    // Derive a channel of normalised matrices
    // - from differential analysis for RNASeq
    // - from normalisedvalidated assays for Affy and MaxQuant
    // - from validated assays for GEO soft file

    ch_norm = ch_differential_norm
        .filter{meta, matrix -> meta.params.study_type == 'rnaseq'}
        .mix(ch_multi_validated_assays.normalised)
        .mix(ch_validated_assays.filter{meta, assay -> meta.params.study_type == 'geo_soft_file'})

    // Prepare channel with normalized matrix, and variance stabilized matrices when available

    ch_processed_matrices = ch_norm.join(ch_differential_varstab, remainder: true)
        .map { meta, norm, vs ->
            def matrices = vs ? [norm, vs].flatten() : [norm]
            [meta, matrices]
        }

    // ========================================================================
    // Functional analysis
    // ========================================================================

    // Prepare background file - for the moment it is only needed for gprofiler2
    // The background file might come from:
    // - the filter matrix, if auto
    // - the gprofiler2 background file, if provided
    // - empty, if not gprofiler2

    ch_background = ch_filtered_matrix
        .filter{meta, matrix -> meta.params.functional_method == 'gprofiler2' && meta.params.gprofiler2_background_file == "auto"}
        .mix(
            ch_paramsets
                .filter{meta -> meta.params.functional_method == 'gprofiler2' && meta.params.gprofiler2_background_file && meta.params.gprofiler2_background_file != "auto"}
                .map{meta -> [meta, file(meta.params.gprofiler2_background_file, checkIfExists: true)]}
        )
        .mix(
            ch_paramsets
                .filter{ meta -> meta.params.functional_method != 'gprofiler2'}
                .map{meta -> [meta, []]}
        )
    // Prepare input for functional analysis

    // - use normalized matrix, if method is gsea
    // - use filtered differential results, if method is gprofiler2
    // - use unfiltered differential results, if method is decoupler
    ch_functional_analysis_matrices = ch_norm
        .filter{meta, matrix -> meta.params.functional_method == 'gsea'}
        .map{ meta, matrix -> [meta, meta, matrix]}
        .mix(
            ch_differential_results_filtered
                // Remember that the meta from differential results contain contrast info.
                // Here the key is the meta without contrast info (same as the meta in other channels)
                // So we can use this key to combine channels
                .filter{meta, meta_with_contrast, results -> meta.params.functional_method == 'gprofiler2'}
        )
        .mix(
            ch_differential_results
                // For decoupler, use unfiltered differential results
                .filter{meta, meta_with_contrast, results -> meta.params.functional_method == 'decoupler'}
        )

    ch_functional_input = ch_functional_analysis_matrices  // meta, meta with contrast, file
        .combine(ch_gene_sets, by:0)             // meta, [gmt files]
        .combine(ch_background, by:0)            // meta, background
        .combine(ch_contrasts, by:0)             // meta, [contrast], [variable], [reference], [target], [formula], [comparison]
        .combine(ch_validated_samplemeta, by:0)  // meta, samplesheet
        .combine(ch_validated_featuremeta, by:0) // meta, features
        .map { it.tail() } // remove the meta without contrast info that is used as key
                        // the remaining meta will contain the contrast info, if any

    // Use a multiMap to generate synched channels for functional enrichment analysis
    ch_functional_input = prepareModuleInput(ch_functional_input, 'functional')
        .multiMap{ meta, input, gene_sets, background, contrasts, variable, reference, target, formula, comparison, samplesheet, features ->
            input: [meta, input, gene_sets, background, meta.params.functional_method]
            contrasts: [meta, contrasts, variable, reference, target, formula, comparison]
            samplesheet: [meta, samplesheet]
            features: [meta, features, meta.params.features_id_col, meta.params.features_name_col]
        }

    // Run functional analysis

    DIFFERENTIAL_FUNCTIONAL_ENRICHMENT(
        ch_functional_input.input,
        ch_functional_input.contrasts,
        ch_functional_input.samplesheet,
        ch_functional_input.features
    )

    // Collect functional analysis results

    ch_functional_results = DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_plot_html
        .join(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_all_enrich, remainder: true)
        .join(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_sub_enrich, remainder: true)
        .mix(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gsea_report)

    // Note that 'functional_method' is additionally added to the meta in the functional subworkflow.
    // Remove it to keep a consistent meta structure (needed for join/combine).
    // Also note that these channels, the meta contain other info apart from the base paramset meta.
    // Hence we create a key using the simple paramset meta with only meta.paramset_name and meta.params,
    // by setting 'use_meta_key' to true. This will facilitate later on to join/combine channels.
    ch_functional_results = prepareModuleOutput(ch_functional_results, ch_paramsets, meta_keys_to_remove=['functional_method'], use_meta_key=true) // key, meta, [ functional results ]

    ch_versions = ch_versions
        .mix(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.versions)


    // ========================================================================
    // Plot figures
    // ========================================================================

    // For geoquery we've done no matrix processing and been supplied with the
    // normalised matrix, which can be passed through to downstream analysis

    ch_mat = ch_norm.filter{meta, matrix -> meta.params.study_type == 'geo_soft_file'}
        .map{meta, norm -> [meta, [norm]]}
        .mix(
            ch_raw
                .filter{meta, raw -> meta.params.study_type != 'geo_soft_file'}
                .join(ch_processed_matrices)
                .map { meta, raw, matrices ->
                    [meta, [raw] + matrices]
                }
        )

    ch_all_matrices = ch_validated_samplemeta                // meta, samples
        .join(ch_validated_featuremeta)                      // meta, features
        .join(ch_mat)                                        // meta, list of matrices (raw, norm, variance stabilized)
        .map { meta, samples, features, matrices ->
            return [meta, samples, features, matrices]
        }

    // Exploratory analysis

    ch_exploratory_input = ch_contrast_variables         // [meta, variable]
        .combine(ch_all_matrices, by:0)                  // [meta, samples, features, [matrices]]
        .map { meta, variable, samples, features, matrices ->
            // we need to add variable info into meta, otherwise the channel will have the same meta no
            // matter the variable information and prepareModuleInput will group them together.
            def meta_new = meta + [id: variable]
            [meta_new, samples, features, matrices, variable]
        }

    PLOT_EXPLORATORY(
        prepareModuleInput(ch_exploratory_input, 'exploratory')
    )

    // Plot differential analysis results

    ch_plot_differential_input = ch_differential_results  // [meta, meta with contrast, results]
        .combine(ch_all_matrices, by: 0)                  // [meta, samples, features, matrices]
        .map{it.tail()}  // remove the meta without contrast info that is used as key
                        // the remaining meta will contain the contrast info

    ch_plot_differential_input = prepareModuleInput(ch_plot_differential_input, 'differential')
        .multiMap{ meta, differential_results, samples, features, matrices ->
            differential_results: [meta, differential_results]
            samples_features_matrices: [meta, samples, features, matrices]
        }

    PLOT_DIFFERENTIAL(
        ch_plot_differential_input.differential_results,
        ch_plot_differential_input.samples_features_matrices
    )

    // Gather software versions

    ch_versions = ch_versions
        .mix(VALIDATOR.out.versions)
        .mix(PLOT_EXPLORATORY.out.versions)
        .mix(PLOT_DIFFERENTIAL.out.versions)

    // ========================================================================
    // ShinyNGS app
    // ========================================================================

    // Make (and optionally deploy) the shinyngs app

    // To prepare the input for shinyngs app we need to first make sure that the differential
    // results are parsed in the same order as the contrast file
    // To do so, as we cannot rely on the order of the channels as they are asynchronous, we
    // create temporary contrast files with the entries based on the order of the gathered
    // differential results

    // As contrasts having 'formula' and 'comparison' won't have a variable,
    // and it is needed for "checkListIsSubset()" in `SHINYNGS` make_app_from_files.R
    // for backwards-compatibility, keep only the channels that have non-empty variable

    // Create a channel with the differential results and the corresponding map with
    // the contrast entries
    differential_with_contrast = ch_paramsets
        .join( ch_differential_results
            .filter { meta, contrast, results -> contrast.variable?.trim() }
            .groupTuple()
        )   // [meta, [meta with contrast], [differential results]]
        .join( ch_contrasts )   // [meta, [contrast], [variable], [reference], [target], [formula], [comparison]]
        .map { meta, meta_with_contrast, results, contrast, variable, reference, target, formula, comparison ->
            // extract the contrast entries from the meta dynamically
            // in this way we don't need to harcode the contrast keys
            def paramset_contrast_keys = contrast[0].keySet()
            def contrast_maps = meta_with_contrast.collect { it.subMap(paramset_contrast_keys) }
            [meta, meta_with_contrast, results, contrast_maps]
        }
        .multiMap { meta, meta_with_contrast, results, contrast_maps ->
            differential_results: [meta, meta_with_contrast, results]
            contrast_maps: [meta, contrast_maps]
        }

    // Save temporary contrast csv files with the entries ordered by the differential results
    ch_contrasts_sorted = differential_with_contrast.contrast_maps
        .collectFile { meta, contrast_map ->
            def header = contrast_map[0].keySet().join(',')
            def content = contrast_map.collect { it.values().join(',') }
            def lines = header + '\n' + content.join('\n') + '\n'
            ["${meta.paramset_name}.csv", lines]
        }
        // parse the channel to have the contrast file with the corresponding meta
        .map { [it.baseName, it] }
        .join( ch_paramsets.map { [it.paramset_name, it] } )
        .map { paramset_name, contrast_file, meta ->
            [meta, contrast_file]
        }

    // Parse input for shinyngs app
    ch_shinyngs_input = differential_with_contrast.differential_results
        .join(ch_contrasts_sorted)
        .join(ch_all_matrices)
        .filter { row ->
            row[0].params.shinyngs_build_app
        }
        .multiMap { meta, meta_with_contrast, differential_results, contrast_file, samplesheet, features, matrices ->
            matrices: [meta, samplesheet, features, matrices]
            contrasts_and_differential: [meta, contrast_file, differential_results]
            contrast_stats_assay: meta.params.exploratory_assay_names.split(',').findIndexOf { it == meta.params.exploratory_final_assay } + 1
        }
    SHINYNGS_APP(
        ch_shinyngs_input.matrices,    // meta, samples, features, [  matrices ]
        ch_shinyngs_input.contrasts_and_differential,   // meta, contrast file, [ differential results ]
        ch_shinyngs_input.contrast_stats_assay
    )

    ch_versions = ch_versions.mix(SHINYNGS_APP.out.versions)

    // ========================================================================
    // Generate report
    // ========================================================================

    // Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'differentialabundance_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'collated_versions.yml', sort: true, newLine: true)

    // Get the report file, logo, css, and citations

    ch_report_files = ch_paramsets
        .map { meta ->
            [ meta, [
                file(meta.params.report_file, checkIfExists: true),
                file(meta.params.logo_file, checkIfExists: true),
                file(meta.params.css_file, checkIfExists: true),
                file(meta.params.citations_file, checkIfExists: true)
            ]]
        }

    // Group differential and functional results by paramset meta [id: study_name, paramset_name: paramset_name, params: paramset]
    // So all the files generated with the same paramset meta will go together to report.
    // Note that the files generated with different contrasts will also be grouped together for the same paramset meta.

    ch_differential_grouped = differential_with_contrast.differential_results.transpose()
        .join(ch_differential_model, by:[0,1])
        .groupTuple()                                 // [ meta, [meta with contrast], [differential results], [differential model] ]
        .map { [it[0], it.tail().tail().flatten()] }  // [ meta, [differential results and models] ]

    ch_functional_grouped = ch_functional_results
        .groupTuple()                                 // [ meta, [meta with contrast], [functional results] ]
        .map { [it[0], it.tail().tail().flatten()] }  // [ meta, [functional results] ]

    // Prepare input for report generation
    // Each paramset will generate one markdown report by gathering all the files created with the same paramset

    ch_report_input = ch_report_files    // [meta, [report_file, logo_file, css_file, citations_file]]
        .combine(ch_collated_versions)   // [versions file]
        .join(ch_all_matrices)           // [meta, samplesheet, features, [matrices]]
        .join(ch_contrasts_sorted)       // [meta, contrast file]
        .join(ch_differential_grouped)   // [meta, [differential results and models]]
        .join(ch_functional_grouped, remainder: true) // [meta, [functional results]]
        .map { [it[0], it.tail().flatten().grep()] }  // [meta, [files]]   // note that grep() would remove null files from join with remainder true
        .map { meta, files -> [meta, files[0], files.tail()] }   // [meta, report_file, [files]]
        .multiMap { meta, report_file, files ->
            report_file:
            [meta, report_file]

            report_params:
            def paramset = [paramset_name: meta.paramset_name] + meta.params  // flat the map
            def report_file_names = ['logo','css','citations','versions_file','observations','features'] +
                paramset.exploratory_assay_names.split(',').collect { "${it}_matrix".toString() } +
                [ 'contrasts_file' ]
            paramset + [report_file_names, files.collect{ f -> f.name}].transpose().collectEntries()

            input_files:
            [meta, files]
        }

    // Render the final report
    if (!params.skip_reports) {
        RMARKDOWNNOTEBOOK(
            ch_report_input.report_file,
            ch_report_input.report_params,
            ch_report_input.input_files.map{ meta, files -> files }
        )

        // Make a report bundle comprising the markdown document and all necessary
        // input files

        ch_bundle_input = RMARKDOWNNOTEBOOK.out.parameterised_notebook
                .combine(ch_report_input.input_files, by:0)
                .map{[it[0], it.tail().flatten()]}   // [meta, [files]]

        MAKE_REPORT_BUNDLE( ch_bundle_input )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

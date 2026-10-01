#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/differentialabundance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/differentialabundance
    Website: https://nf-co.re/differentialabundance
    Slack  : https://nfcore.slack.com/channels/differentialabundance
----------------------------------------------------------------------------------------
*/

nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIALABUNDANCE   } from './workflows/differentialabundance'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'
include { buildParamset           } from './subworkflows/local/utils_nfcore_differentialabundance_pipeline'

// A file or set of files that the pipeline publishes, with the name that decides where it goes
record Published {
    name:  String
    meta:  Map
    files: List<Path>
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The files of a run (input, contrasts, matrix, feature_length_matrix, gtf) are optional values, so a
    pipeline that includes this one can give them from its own dataflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Input options
    input: Value<Path>?

    study_name: String = 'study'
    study_type: String = 'rnaseq'
    study_abundance_type: String = 'abundance'
    contrasts: Value<Path>?
    querygse: String?
    matrix: Value<Path>?
    feature_length_matrix: Value<Path>?
    control_features: String?
    sizefactors_from_controls: Boolean = false

    // Output options
    round_digits: Integer = -1
    seed: Integer?

    // Paramsheet options
    paramset_name: String?
    paramsheet: String?

    // Sample sheet options
    observations_type: String = 'sample'
    observations_id_col: String = 'sample'
    observations_name_col: String?

    // Feature options
    features: String?
    features_type: String = 'gene'
    features_id_col: String = 'gene_id'
    features_name_col: String = 'gene_name'
    features_metadata_cols: String = 'gene_id,gene_name,gene_biotype'

    // GTF parsing options
    gtf: Value<Path>?
    features_gtf_feature_type: String = 'transcript'
    features_gtf_table_first_field: String = 'gene_id'

    // Affy-specific options
    affy_cel_files_archive: String?
    affy_file_name_col: String = 'file'
    affy_background: Boolean = true
    affy_bgversion: Integer = 2
    affy_destructive: Boolean = false
    affy_cdfname: String?
    affy_rm_mask: Boolean = false
    affy_rm_outliers: Boolean = false
    affy_rm_extra: Boolean = false
    affy_build_annotation: Boolean = true

    // Proteus-specific options
    proteus_measurecol_prefix: String = 'LFQ intensity'
    proteus_norm_function: String = 'normalizeMedian'
    proteus_plotsd_method: String = 'violin'
    proteus_plotmv_loess: Boolean = true
    proteus_palette_name: String = 'Set1'

    // Filtering options
    filtering_min_samples: Float = 1.0
    filtering_min_abundance: Object = 1
    filtering_min_proportion: Float?
    filtering_grouping_var: String?
    filtering_min_proportion_not_na: Float = 0.5
    filtering_min_samples_not_na: Float?

    // Exploratory options
    exploratory_main_variable: String = 'auto_pca'
    exploratory_clustering_method: String = "ward.D2"
    exploratory_cor_method: String = "spearman"
    exploratory_n_features: Integer = 500
    exploratory_whisker_distance: Float = 1.5
    exploratory_mad_threshold: Integer = -5
    exploratory_assay_names: String = "raw,normalised,variance_stabilised"
    exploratory_final_assay: String = "variance_stabilised"
    exploratory_log2_assays: String? = 'raw,normalised'
    exploratory_palette_name: String = 'Set1'

    // Differential options
    differential_method: String = 'deseq2'  // 'deseq2', 'limma', 'dream', 'propd'
    differential_file_suffix: String?
    differential_feature_id_column: String = "gene_id"
    differential_min_fold_change: Float = 2.0
    differential_max_pval: Float = 1.0
    differential_max_qval: Float = 0.05
    differential_palette_name: String = 'Set1'
    differential_subset_to_contrast_samples: Boolean = false

    // DESeq2-specific options
    deseq2_test: String = "Wald"
    deseq2_fit_type: String = "parametric"
    deseq2_sf_type: String = 'ratio'
    deseq2_min_replicates_for_replace: Integer = 7
    deseq2_use_t: Boolean = false
    deseq2_lfc_threshold: Float = 0
    deseq2_alt_hypothesis: String = 'greaterAbs'
    deseq2_independent_filtering: Boolean = true
    deseq2_p_adjust_method: String = 'BH'
    deseq2_alpha: Float = 0.1
    deseq2_minmu: Float = 0.5
    deseq2_vs_method: String = 'vst'  // 'rlog', 'vst', or 'rlog,vst'
    deseq2_shrink_lfc: Boolean = true
    deseq2_vs_blind: Boolean = true
    deseq2_vst_nsub: Integer = 1000

    // Limma-specific options
    limma_ndups: Float?
    limma_spacing: String?
    limma_block: String?
    limma_correlation: String?
    limma_method: String = 'ls'
    limma_proportion: Float = 0.01
    limma_stdev_coef_lim: String = '0.1,4'
    limma_trend: Boolean = false
    limma_robust: Boolean = false
    limma_winsor_tail_p: String = '0.05,0.1'
    limma_adjust_method: String = "BH"
    limma_p_value: Float = 1.0
    limma_lfc: Integer = 0
    limma_confint: Boolean = false
    limma_use_voom: Boolean = false

    // DREAM-specific options
    dream_adjust_method: String = "BH"
    dream_p_value: Integer = 1
    dream_lfc: Integer = 0
    dream_confint: Boolean = false
    dream_proportion: Float = 0.01
    dream_stdev_coef_lim: String = '0.1,4'
    dream_trend: Boolean = false
    dream_robust: Boolean = false
    dream_winsor_tail_p: String = '0.05,0.1'
    dream_ddf: String = 'adaptive'
    dream_reml: Boolean = false
    dream_apply_voom: Boolean = false

    // propd-specific options
    propd_alpha: Float?
    propd_moderated: Boolean = true
    propd_fdr: Float = 0.05
    propd_permutation: Integer = 0
    propd_number_of_cutoffs: Integer = 100
    propd_save_pairwise: Boolean = false
    propd_save_pairwise_full: Boolean = false
    propd_save_adjacency: Boolean = false
    propd_save_rdata: Boolean = false

    // functional analysis options
    functional_method: String = 'none'  // 'none', 'gsea', 'gprofiler2', 'decoupler', 'grea'
    gene_sets_files: String?

    // GSEA options
    gsea_nperm: Integer = 1000
    gsea_permute: String = 'phenotype'
    gsea_scoring_scheme: String = 'weighted'
    gsea_metric: String = 'Signal2Noise'
    gsea_sort: String = 'real'
    gsea_order: String = 'descending'
    gsea_set_max: Integer = 500
    gsea_set_min: Integer = 15

    gsea_norm: String = 'meandiv'
    gsea_rnd_type: String = 'no_balance'
    gsea_make_sets: Boolean = true
    gsea_median: Boolean = false
    gsea_num: Integer = 100
    gsea_plot_top_x: Integer = 20
    gsea_save_rnd_lists: Boolean = false
    gsea_zip_report: Boolean = false

    // gprofiler2 options
    gprofiler2_organism: String?
    gprofiler2_significant: Boolean = true
    gprofiler2_measure_underrepresentation: Boolean = false
    gprofiler2_correction_method: String = 'gSCS'
    gprofiler2_sources: String?
    gprofiler2_evcodes: Boolean = false
    gprofiler2_max_qval: Float = 0.05
    gprofiler2_token: String?
    gprofiler2_background_file: String = 'auto'
    gprofiler2_background_column: String?
    gprofiler2_domain_scope: String = 'annotated'
    gprofiler2_min_diff: Integer = 1
    gprofiler2_palette_name: String = 'Blues'

    // decoupler options
    decoupler_network: String?
    decoupler_min_n: Integer = 5
    decoupler_methods: String = 'ulm'

    // grea options
    grea_set_min: Integer = 15
    grea_set_max: Integer = 500
    grea_permutation: Integer = 100

    // ShinyNGS
    shinyngs_build_app: Boolean = true

    // Report options
    skip_reports: Boolean = false

    // Note: for shinyapps deployment, in addition to setting these values,
    // SHINYAPPS_TOKEN and SHINYAPPS_SECRET must be available to the
    // environment, probably via Nextflow secrets
    shinyngs_deploy_to_shinyapps_io: Boolean = false
    shinyngs_shinyapps_account: String?
    shinyngs_shinyapps_app_name: String?

    // Reporting
    logo_file: String = "${moduleDir}/docs/images/nf-core-differentialabundance_logo_light.png"
    css_file: String = "${moduleDir}/assets/nf-core_style.css"
    citations_file: String = "${moduleDir}/CITATIONS.md"
    report_file: String = "${moduleDir}/assets/differentialabundance_report.qmd"
    report_title: String?
    report_author: String?
    report_contributors: String?
    report_description: String?
    disable_report_modules: String?

    // References
    genome: String?

    // Boilerplate options
    email: String?
    email_on_fail: String?
    plaintext_email: Boolean = false
    monochrome_logs: Boolean = false
    help: Boolean = false
    help_full: Boolean = false
    show_hidden: Boolean = false
    version: Boolean = false

    // Schema validation default options
    validate_params: Boolean = true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks. The paramsets of a paramsheet come from here.
    //
    init = PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.help,
        params.help_full,
        params.show_hidden,
        params
    )

    //
    // Without a paramsheet there is one paramset, built from the params and from the files of the run
    //
    def ch_null = channel.value(null as Path)
    def genome_gtf = getGenomeAttribute('gtf')
    def ch_gtf: Value<Path> = params.gtf != null ? params.gtf : genome_gtf != null ? channel.value(file(genome_gtf)) : ch_null

    def ch_files = (params.input != null ? params.input : ch_null)
        .map { input -> record(input: input) }
        .combine(contrasts: params.contrasts != null ? params.contrasts : ch_null)
        .combine(matrix: params.matrix != null ? params.matrix : ch_null)
        .combine(feature_length_matrix: params.feature_length_matrix != null ? params.feature_length_matrix : ch_null)
        .combine(gtf: ch_gtf)

    def ch_paramsets = init.paramsets
    if (!params.paramsheet) {
        ch_paramsets = ch_files.flatMap { files -> [ buildParamset(params, files) ] }
    }

    //
    // WORKFLOW: Run main workflow
    //
    abundance = DIFFERENTIALABUNDANCE (
        ch_paramsets
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.monochrome_logs
    )

    //
    // Build category channels for publishing
    //
    ch_pub_preprocessing = abundance.affy_cel_files.map { item -> record(name: 'affy_cel_files', meta: item[0], files: item[1..-1]) }
        .mix(abundance.affy_raw_expression.map        { meta, file -> record(name: 'affy_raw_expression', meta: meta, files: [file].flatten()) })
        .mix(abundance.affy_norm_expression.map       { meta, file -> record(name: 'affy_norm_expression', meta: meta, files: [file].flatten()) })
        .mix(abundance.affy_annotation.map            { meta, file -> record(name: 'affy_annotation', meta: meta, files: [file].flatten()) })
        .mix(abundance.affy_raw_rds.map               { meta, file -> record(name: 'affy_raw_rds', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_raw.map                { meta, file -> record(name: 'proteus_raw', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_norm.map               { meta, file -> record(name: 'proteus_norm', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_plots.map              { meta, file -> record(name: 'proteus_plots', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_raw_rdata.map          { meta, file -> record(name: 'proteus_raw_rdata', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_norm_rdata.map         { meta, file -> record(name: 'proteus_norm_rdata', meta: meta, files: [file].flatten()) })
        .mix(abundance.proteus_session_info.map       { meta, file -> record(name: 'proteus_session_info', meta: meta, files: [file].flatten()) })
        .mix(abundance.geo_expression.map             { meta, file -> record(name: 'geo_expression', meta: meta, files: [file].flatten()) })
        .mix(abundance.geo_annotation.map             { meta, file -> record(name: 'geo_annotation', meta: meta, files: [file].flatten()) })
        .mix(abundance.geo_rds.map                    { meta, file -> record(name: 'geo_rds', meta: meta, files: [file].flatten()) })
        .mix(abundance.gtf_annotation.map             { meta, file -> record(name: 'gtf_annotation', meta: meta, files: [file].flatten()) })

    ch_pub_differential = abundance.diff_results.map              { _key, meta, file -> record(name: 'results', meta: meta, files: [file].flatten()) }
        .mix(abundance.diff_results_filtered.map  { _key, meta, file -> record(name: 'results_filtered', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_normalised_matrix.map { meta, file -> record(name: 'normalised_matrix', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_variance_stabilised.map { meta, file -> record(name: 'variance_stabilised_matrix', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_size_factors.map      { meta, file -> record(name: 'size_factors', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_dispersion_plot.map   { meta, file -> record(name: 'dispersion_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_md_plot.map           { meta, file -> record(name: 'md_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_rdata.map             { meta, file -> record(name: 'rdata', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_session_info.map      { meta, file -> record(name: 'session_info', meta: meta, files: [file].flatten()) })
        .mix(abundance.diff_annotated.map         { meta, file -> record(name: 'annotated', meta: meta, files: [file].flatten()) })

    ch_pub_functional = abundance.gsea_report_tsv.map          { meta, ref, target -> record(name: 'gsea_report_tsv', meta: meta, files: [ref, target].flatten()) }
        .mix(abundance.gsea_report_html.map       { meta, ref, target -> record(name: 'gsea_report_html', meta: meta, files: [ref, target].flatten()) })
        .mix(abundance.gsea_index_html.map        { meta, file -> record(name: 'gsea_index_html', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_heat_map_corr_plot.map { meta, file -> record(name: 'gsea_heat_map_corr_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_ranked_gene_list.map  { meta, file -> record(name: 'gsea_ranked_gene_list', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_sizes.map    { meta, file -> record(name: 'gsea_gene_set_sizes', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_histogram.map         { meta, file -> record(name: 'gsea_histogram', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_heatmap.map           { meta, file -> record(name: 'gsea_heatmap', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_pvalues_vs_nes_plot.map { meta, file -> record(name: 'gsea_pvalues_vs_nes_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_ranked_list_corr.map  { meta, file -> record(name: 'gsea_ranked_list_corr', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_butterfly_plot.map    { meta, file -> record(name: 'gsea_butterfly_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_tsv.map      { meta, file -> record(name: 'gsea_gene_set_tsv', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_html.map     { meta, file -> record(name: 'gsea_gene_set_html', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_heatmap.map  { meta, file -> record(name: 'gsea_gene_set_heatmap', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_enplot.map   { meta, file -> record(name: 'gsea_gene_set_enplot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_gene_set_dist.map     { meta, file -> record(name: 'gsea_gene_set_dist', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_snapshot.map          { meta, file -> record(name: 'gsea_snapshot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_archive.map           { meta, file -> record(name: 'gsea_archive', meta: meta, files: [file].flatten()) })
        .mix(abundance.gsea_rpt.map               { meta, file -> record(name: 'gsea_rpt', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_html.map        { meta, file -> record(name: 'gprofiler2_html', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_all_enrichment.map { meta, file -> record(name: 'gprofiler2_all_enrichment', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_sub_enrichment.map { meta, file -> record(name: 'gprofiler2_sub_enrichment', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_plot_png.map    { meta, file -> record(name: 'gprofiler2_plot_png', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_sub_plot.map    { meta, file -> record(name: 'gprofiler2_sub_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_rds.map         { meta, file -> record(name: 'gprofiler2_rds', meta: meta, files: [file].flatten()) })
        .mix(abundance.gprofiler2_filtered_gmt.map { meta, file -> record(name: 'gprofiler2_filtered_gmt', meta: meta, files: [file].flatten()) })
        .mix(abundance.decoupler_estimate.map     { meta, file -> record(name: 'decoupler_estimate', meta: meta, files: [file].flatten()) })
        .mix(abundance.decoupler_pvals.map        { meta, file -> record(name: 'decoupler_pvals', meta: meta, files: [file].flatten()) })
        .mix(abundance.decoupler_png.map          { meta, file -> record(name: 'decoupler_png', meta: meta, files: [file].flatten()) })
        .mix(abundance.functional_session_info.map { meta, file -> record(name: 'session_info', meta: meta, files: [file].flatten()) })

    ch_pub_plotting = abundance.plot_exploratory.map { meta, file -> record(name: 'exploratory', meta: meta, files: [file].flatten()) }
        .mix(abundance.plot_volcanos.map          { meta, file -> record(name: 'differential_volcanos', meta: meta, files: [file].flatten()) })

    ch_pub_shinyngs = abundance.shinyngs_data.map { meta, file -> record(name: 'shinyngs_data', meta: meta, files: [file].flatten()) }
        .mix(abundance.shinyngs_app_file.map      { meta, file -> record(name: 'shinyngs_app', meta: meta, files: [file].flatten()) })

    ch_pub_report = abundance.report_html.map { meta, file -> record(name: 'report_html', meta: meta, files: [file].flatten()) }
        .mix(abundance.report_bundle.map          { meta, file -> record(name: 'report_bundle', meta: meta, files: [file].flatten()) })

    ch_pub_versions = abundance.nfcore_versions.map { file -> record(name: 'versions', meta: [:], files: [file].flatten()) }
        .mix(abundance.collated_versions.map      { file -> record(name: 'collated_versions', meta: [:], files: [file].flatten()) })

    // The files of propd and grea go in folders named after the run, without the paramset
    ch_pub_propd = abundance.propd_results.map  { meta, file -> record(name: 'results', meta: meta, files: [file].flatten()) }
        .mix(abundance.propd_pairwise.map          { meta, file -> record(name: 'pairwise', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_pairwise_filtered.map { meta, file -> record(name: 'pairwise_filtered', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_fdr.map               { meta, file -> record(name: 'fdr', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_genewise_plot.map     { meta, file -> record(name: 'genewise_plot', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_rdata.map             { meta, file -> record(name: 'rdata', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_adjacency.map         { meta, file -> record(name: 'adjacency', meta: meta, files: [file].flatten()) })
        .mix(abundance.propd_session_info.map      { file -> record(name: 'session_info', meta: [id: file.name.replace('.R_sessionInfo.log', '')], files: [file].flatten()) })

    ch_pub_grea = abundance.grea_results.map { meta, file -> record(name: 'results', meta: meta, files: [file].flatten()) }
        .mix(abundance.grea_session_info.map { file -> record(name: 'session_info', meta: [id: file.name.replace('.R_sessionInfo.log', '')], files: [file].flatten()) })

    publish:
    preprocessing = ch_pub_preprocessing
    differential  = ch_pub_differential
    functional    = ch_pub_functional
    propd         = ch_pub_propd
    grea          = ch_pub_grea
    plotting      = ch_pub_plotting
    shinyngs      = ch_pub_shinyngs
    report        = ch_pub_report
    versions      = ch_pub_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISHING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def preprocessingTarget(r: Published) -> String {
    def folder = [
        // AFFY
        affy_raw_expression        : 'tables/processed_abundance',
        affy_norm_expression       : 'tables/processed_abundance',
        affy_annotation            : 'tables/annotation',
        affy_raw_rds               : 'other/affy',
        affy_cel_files             : 'untar',

        // PROTEUS
        proteus_raw                : 'tables/proteus',
        proteus_norm               : 'tables/proteus',
        proteus_plots              : 'plots/proteus',
        proteus_raw_rdata          : 'other/proteus',
        proteus_norm_rdata         : 'other/proteus',
        proteus_session_info       : 'other/proteus',

        // GEO SOFT
        geo_expression             : 'tables/processed_abundance',
        geo_annotation             : 'tables/annotation',
        geo_rds                    : 'other/affy',

        // GTF
        gtf_annotation             : 'tables/annotation',

    ][r.name] ?: r.name
    def target = (r.name in ['proteus_plots', 'proteus_raw_rdata', 'proteus_norm_rdata']) \
        ? "${folder}/${r.meta.paramset_name}/${r.meta.contrast}/" \
        : "${folder}/${r.meta.paramset_name}/"
    return target
}

def differentialTarget(r: Published) -> String {
    def folder = [
        results                    : 'tables/differential',
        results_filtered           : 'tables/differential',
        annotated                  : 'tables/differential',
        normalised_matrix          : 'tables/processed_abundance',
        variance_stabilised_matrix : 'tables/processed_abundance',
        size_factors               : "other/${r.meta.params.differential_method}",
        dispersion_plot            : 'plots/qc',
        md_plot                    : 'plots/qc',
        rdata                      : "other/${r.meta.params.differential_method}",
        session_info               : "other/${r.meta.params.differential_method}"
    ][r.name] ?: r.name
    return "${folder}/${r.meta.paramset_name}/"
}

def propdTarget(r: Published) -> String {
    def folder = [
        results           : 'tables/differential',
        pairwise          : 'tables/differential',
        pairwise_filtered : 'tables/differential',
        fdr               : 'tables/differential',
        genewise_plot     : "plots/differential/${r.meta.id}",
        rdata             : "other/propd/${r.meta.id}",
        adjacency         : "other/propd/${r.meta.id}",
        session_info      : "other/propd/${r.meta.id}"
    ][r.name] ?: r.name
    return "${folder}/"
}

def greaTarget(r: Published) -> String {
    def folder = [
        results      : "tables/functional/grea/${r.meta.id}",
        session_info : "other/grea/${r.meta.id}"
    ][r.name] ?: r.name
    return "${folder}/"
}

def functionalTarget(r: Published) -> String {
    def method = r.meta.params.functional_method
    def folder = [
        // GSEA
        gsea_report_tsv           : 'report/gsea',
        gsea_report_html          : 'report/gsea',
        gsea_index_html           : 'report/gsea',
        gsea_heat_map_corr_plot   : 'report/gsea',
        gsea_ranked_gene_list     : 'report/gsea',
        gsea_gene_set_sizes       : 'report/gsea',
        gsea_histogram            : 'report/gsea',
        gsea_heatmap              : 'report/gsea',
        gsea_pvalues_vs_nes_plot  : 'report/gsea',
        gsea_ranked_list_corr     : 'report/gsea',
        gsea_butterfly_plot       : 'report/gsea',
        gsea_gene_set_tsv         : 'report/gsea',
        gsea_gene_set_html        : 'report/gsea',
        gsea_gene_set_heatmap     : 'report/gsea',
        gsea_gene_set_enplot      : 'report/gsea',
        gsea_gene_set_dist        : 'report/gsea',
        gsea_snapshot             : 'report/gsea',
        gsea_archive              : 'report/gsea',
        gsea_rpt                  : 'report/gsea',

        // GPROFILER2
        gprofiler2_all_enrichment : 'tables/gprofiler2',
        gprofiler2_sub_enrichment : 'tables/gprofiler2',
        gprofiler2_html           : 'plots/gprofiler2',
        gprofiler2_plot_png       : 'plots/gprofiler2',
        gprofiler2_sub_plot       : 'plots/gprofiler2',
        gprofiler2_rds            : 'other/gprofiler2',
        gprofiler2_filtered_gmt   : 'other/gprofiler2',

        // DECOUPLER
        decoupler_estimate        : 'tables/decoupler',
        decoupler_pvals           : 'tables/decoupler',
        decoupler_png             : 'plots/decoupler',

        // common outputs
        session_info               : "other/${method}"
    ][r.name] ?: r.name
    def gene_set_name = (method == 'gsea' && r.meta.params.gene_sets_files) \
        ? r.meta.params.gene_sets_files.tokenize('/')[-1].replaceFirst(/\.[^.]+$/, '') \
        : null
    def target = (method == 'gsea') ? "${folder}/${r.meta.paramset_name}/${r.meta.id}/${gene_set_name}/" \
        : (method == 'gprofiler2') ? "${folder}/${r.meta.paramset_name}/${r.meta.id}/" \
        : "${folder}/${r.meta.paramset_name}/"
    return target
}

def plottingTarget(r: Published) -> String {
    def folder = [
        exploratory           : 'plots/exploratory',
        differential_volcanos : 'plots/differential',
    ][r.name] ?: r.name
    return "${folder}/${r.meta.paramset_name}/"
}

output {
    preprocessing: Channel<Published> {
        path { r -> r.files >> preprocessingTarget(r) }
    }
    differential: Channel<Published> {
        path { r -> r.files >> differentialTarget(r) }
    }
    functional: Channel<Published> {
        path { r -> r.files >> functionalTarget(r) }
    }
    plotting: Channel<Published> {
        path { r -> r.files >> plottingTarget(r) }
    }
    propd: Channel<Published> {
        path { r -> r.files >> propdTarget(r) }
    }
    grea: Channel<Published> {
        path { r -> r.files >> greaTarget(r) }
    }
    shinyngs: Channel<Published> {
        path { r -> r.files >> "shinyngs_app/${r.meta.paramset_name}/" }
    }
    report: Channel<Published> {
        path { r -> r.files >> "report/${r.meta.paramset_name}/" }
    }
    versions: Channel<Published> {
        path { r -> r.files >> "pipeline_info/" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

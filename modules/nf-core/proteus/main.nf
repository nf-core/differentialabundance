process PROTEUS {
    tag "$meta"
    label 'process_medium'
//TODO: Change containers
//    conda "proteus2"
//    conda "bioconda::r-proteus-bartongroup=0.2.16 conda-forge::r-plotly=4.10.1 bioconda::bioconductor-limma=3.54.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-proteus-bartongroup:0.2.16--r42hdfd78af_0' :
        'quay.io/biocontainers/mulled-v2-0ad0abd3e3e02e24e1626edaef6d6f9a967733fb:246a3b59c610a9cd35bdf8110be3a908e3769ae0-0' }"

    input:
    tuple val(meta), path(samplesheet), path(quants)

    output:
    tuple val(meta), path("*normalised_distributions.png")      , emit: nonnorm_dist_plot
    tuple val(meta), path("*normalised_distributions.png")      , emit: norm_dist_plot
    tuple val(meta), path("*mean_variance_relationship.png")    , emit: mean_var_relationship_plot
    tuple val(meta), path("*dendrogram.png")                    , emit: dendro_plot
    tuple val(meta), path("*raw_proteingroups.rds")                 , emit: rdata
    tuple val(meta), path("*raw_proteingroups_tab.tsv")             , emit: tab
    tuple val(meta), path("*normalised_proteingroups_tab.tsv")  , emit: normtab
//    tuple val(meta), path("*normalised_proteingroups_tab2.tsv")  , emit: normtab2
    tuple val(meta), path("*R_sessionInfo.log")                , emit: session_info
//    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'proteus_readproteingroups.R'
}

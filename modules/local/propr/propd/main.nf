process PROPR_PROPD {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:75e6dfd62a2313a9':
        'community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:fbd569ab00953cb0' }"

    input:
    tuple val(meta), path(count), path(samplesheet), val(contrast_variable), val(reference), val(target)

    output:
    tuple val(meta), path("*.propd.rds")                 , emit: rds
    tuple val(meta), path("*.propd.results.tsv")         , emit: results
    tuple val(meta), path("*.propd.results_filtered.tsv"), emit: results_filtered, optional: true
    tuple val(meta), path("*.propd.adjacency.csv")       , emit: adjacency       , optional: true
    tuple val(meta), path("*.propd.connectivity.tsv")    , emit: connectivity    , optional: true
    tuple val(meta), path("*.propd.hub_genes.tsv")       , emit: hub_genes       , optional: true
    tuple val(meta), path("*.propd.fdr.tsv")             , emit: fdr             , optional: true
    path "*.warnings.log"                                , emit: warnings
    path "*.R_sessionInfo.log"                           , emit: session_info
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propd.R'
}

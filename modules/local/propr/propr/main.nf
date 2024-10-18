process PROPR_PROPR {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:75e6dfd62a2313a9':
        'community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:fbd569ab00953cb0' }"

    input:
    tuple val(meta), path(count)

    output:
    tuple val(meta), path("*.propr.rds")          , emit: propr
    tuple val(meta), path("*.propr.matrix.csv")   , emit: matrix
    tuple val(meta), path("*.propr.fdr.tsv")      , emit: fdr         , optional:true
    tuple val(meta), path("*.propr.adjacency.csv"), emit: adjacency   , optional:true
    path "*.warnings.log"                         , emit: warnings
    path "*.R_sessionInfo.log"                    , emit: session_info
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propr.R'
}

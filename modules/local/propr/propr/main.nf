process PROPR_PROPR {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-401a215d4024df776a98d90a352048199e342a3d:5ba9bbf6cd4f4f98983526673c223d2e7d829b36-0':
    //     'biocontainers/mulled-v2-401a215d4024df776a98d90a352048199e342a3d:5ba9bbf6cd4f4f98983526673c223d2e7d829b36-0' }"

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

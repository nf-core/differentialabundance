process PROPR_GREA {
    tag "$meta.id"
    label 'process_high'

    // TODO
    // the container is not updated for the grea function
    // for the moment I'm running it locally
    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-401a215d4024df776a98d90a352048199e342a3d:5ba9bbf6cd4f4f98983526673c223d2e7d829b36-0':
    //     'biocontainers/mulled-v2-401a215d4024df776a98d90a352048199e342a3d:5ba9bbf6cd4f4f98983526673c223d2e7d829b36-0' }"

    input:
    tuple val(meta), path(adj)
    tuple val(meta2), path(gmt)

    output:
    tuple val(meta), path("*.grea.tsv"), emit: results
    path "versions.yml",                 emit: versions
    path "*.R_sessionInfo.log",          emit: session_info

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'grea.R'
}

process FILTERVAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-propr:5.0.3':
        'quay.io/biocontainers/r-propr:5.0.3' }"

    input:
    tuple val(meta), path(count), path(adj_matrix)

    output:
    tuple val(meta), path("*.count_filtered.tsv"),  emit: count
    path "*.R_sessionInfo.log",                     emit: session_info
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'filtervar.R'
}
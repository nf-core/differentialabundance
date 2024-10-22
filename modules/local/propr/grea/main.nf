process PROPR_GREA {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:75e6dfd62a2313a9':
        'community.wave.seqera.io/library/bioconductor-limma_r-ggplot2_r-propr:fbd569ab00953cb0' }"

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

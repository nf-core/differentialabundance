process PROPR_PROPR {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/suzannejin/propr:latest"

    input:
    tuple val(meta), path(count)

    output:
    tuple val(meta), path("*.propr.rds"), emit: propr
    tuple val(meta), path("*.propr.tsv"), emit: matrix
    tuple val(meta), path("*.fdr.tsv"),   emit: fdr         , optional:true
    path "*.R_sessionInfo.log",           emit: session_info
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'propr.R'
}

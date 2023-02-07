process IMPORT_FEATURECOUNTS {
    tag "$meta"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::r-base=4.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base%3A4.2.1':
        'quay.io/biocontainers/r-base:4.2.1' }"

    input:
    path(counts)
    path(samplesheet)

    output:
    path("*.tsv")                           , emit: counts
    path("*.csv")                           , emit: samplesheet
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'import_featurecounts.R'
}

process QUARTODOCUMENT {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/giusmar/quarto:1.3.433'

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    /*
    tuple val(meta), path("*.html")           , emit: report
    tuple val(meta), path ("artifacts/*")     , emit: artifacts, optional: true
    tuple val(meta), path ("session_info.log"), emit: session_info
    path  "versions.yml"                      , emit: versions
    */
    path "*"                                  , includeInputs: true

    script:
    """
    echo ciao
    """
}
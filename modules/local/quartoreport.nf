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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parametrize = (task.ext.parametrize == null) ?  true : task.ext.parametrize
    def implicit_params = (task.ext.implicit_params == null) ? true : task.ext.implicit_params
    def meta_params = (task.ext.meta_params == null) ? true : task.ext.meta_params
    """
    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}
    """
}
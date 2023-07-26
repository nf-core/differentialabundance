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

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (parametrize) {
        nb_params = [:]
        if (implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = """\
            params = yaml::read_yaml('.params.yml')
            rmarkdown::render('${prefix}.Rmd', params=params, envir=new.env())
        """
    } else {
        render_cmd = "rmarkdown::render('${prefix}.Rmd')"
    }
    """
    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}
    """
}
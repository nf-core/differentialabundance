process EXECUTE_REPORT {
    tag "report"
    label 'process_medium'

    //    TODO: See how the container is best solved (i.e. provide container name as param or create one huge container containing all possible packages that the different report could require)
//    conda (params.enable_conda ? "conda-forge::r-base=4.1.2 bioconda::r-alakazam=1.2.0 bioconda::r-shazam=1.1.0 conda-forge::r-kableextra=1.3.4 conda-forge::r-knitr=1.33 conda-forge::r-stringr=1.4.0 conda-forge::r-dplyr=1.0.6 conda-forge::r-optparse=1.7.1" : null)
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/mulled-v2-7da73314bcc47157b442d16c3dcfbe81e75a404f:9bb35f8114dffcd97b3afb5de8587355aca16b66-0' :
//        'quay.io/biocontainers/mulled-v2-7da73314bcc47157b442d16c3dcfbe81e75a404f:9bb35f8114dffcd97b3afb5de8587355aca16b66-0' }"

    container 'qbicpipelines/rnadeseq:dev'

    input:
    path(report_file)
    val(args)

    output:
//    path "versions.yml" , emit: versions
    path("*.html"), emit: report_html
    path("*.txt"), emit: txt
    path("*.png"), emit: png
    path("*.pdf"), emit: pdf

//    path(repertoire_report)

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'execute_report.R'
}

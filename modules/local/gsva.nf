process COMPUTE_GSVA{
    label 'process_single'

    container 'biocontainers/mulled-v2-dfd42984f8ba1b1d0d13648f30430581dff51e82:2a27d1707e91426237b882a7ba5adec8724f9069-0'

    input:
    path counts

    output:
    path "*.txt", emit: gsva_for_mqc
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    compute_gsva.R $params.genesets $counts $params.species

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-GSVA: \$(Rscript -e "library(GSVA); cat(as.character(packageVersion('GSVA')))")
        r-stringr: \$(Rscript -e "library(stringr); cat(as.character(packageVersion('stringr')))")
    END_VERSIONS
    """
}

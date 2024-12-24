process VALIDATE_MODEL {
    tag "${models_yml}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f5/f5ef116da0b04fc6359efa48f08b9acfc040a02d7a1ec62476c77b6db0fc1499/data' :
        'community.wave.seqera.io/library/r-jsonlite_r-optparse_r-tidyverse_r-yaml:18dc3fc2d990206d' }"

    input:
    tuple val(meta) , path(samplesheet)
    tuple val(meta2), path(models_yml)

    output:
    tuple val(meta), path("pheno_table.csv"), emit: pheno_table
    tuple val(meta), path("designs.json")   , emit: designs
    tuple val(meta), path("warnings.json")  , emit: warnings, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    validate_model.R \\
        --yml "${models_yml}" \\
        --samplesheet "${samplesheet}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'R.Version()\$version.string' | sed -n 's|\\[1\\] "R version \\(.*\\) (.*|\\1|p')
        tidyverse: \$(Rscript -e "cat(paste(packageVersion('tidyverse'), collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse'), collapse='.'))")
        yaml: \$(Rscript -e "cat(paste(packageVersion('yaml'), collapse='.'))")
        jsonlite: \$(Rscript -e "cat(paste(packageVersion('jsonlite'), collapse='.'))")
    END_VERSIONS
    """

    stub:
    """
    touch pheno_table.csv
    touch designs.json
    touch warnings.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'R.Version()\$version.string' | sed -n 's|\\[1\\] "R version \\(.*\\) (.*|\\1|p')
        tidyverse: \$(Rscript -e "cat(paste(packageVersion('tidyverse'), collapse='.'))")
        optparse: \$(Rscript -e "cat(paste(packageVersion('optparse'), collapse='.'))")
        yaml: \$(Rscript -e "cat(paste(packageVersion('yaml'), collapse='.'))")
        jsonlite: \$(Rscript -e "cat(paste(packageVersion('jsonlite'), collapse='.'))")
    END_VERSIONS
    """

}

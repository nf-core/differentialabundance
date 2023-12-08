process FILTER_DIFFTABLE {

    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(tsv)
    tuple val(logFC_column), val(logFC_threshold)
    tuple val(padj_column), val(padj_threshold)

    output:
    tuple val(meta), path("*_filtered.tsv") , emit: filtered
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '9.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    function find_column_number {
        file=\$1
        column=\$2

        head -n 1 \$file | tr '\\t' '\\n' | grep -n "^\${column}\$" | awk -F':' '{print \$1}'
    }

    # export the two threshold as otherwise, awk will for some reason interpret them incorrectly; padj_threshold 0.05 then becomes 0
    export logFC_threshold
    export padj_threshold


    logFC_threshold=$logFC_threshold
    padj_threshold=$padj_threshold

    logFC_col=\$(find_column_number $tsv $logFC_column)
    padj_col=\$(find_column_number $tsv $padj_column)
    outfile=\$(echo $tsv | sed 's/\\(.*\\)\\..*/\\1/')_filtered.tsv

    head -n 1 $tsv > \${outfile}.tmp

    # The following snippet performs the following checks on each row (add +0.0 to the numbers so that they are definitely treated as numerics):
    # 1. Check that the current logFC/padj is not NA
    # 2. Check that the current logFC is >= threshold (abs does not work, so use a workaround)
    # 3. Check that the current padj is <= threshold
    # If this is true, the row is written to the new file, otherwise not
    tail -n +2 $tsv | awk -F'\\t' -v logFC_col=\$logFC_col -v padj_col=\$padj_col '
    \$logFC_col != "NA" && \$padj_col != "NA" && (\$logFC_col+0.0 >= 0 ? \$logFC_col+0.0 >= ENVIRON["logFC_threshold"]+0.0 : -(\$logFC_col+0.0) >= ENVIRON["logFC_threshold"]+0.0) && \$padj_col+0.0 <= ENVIRON["padj_threshold"]+0.0 {print}' >> \${outfile}.tmp

    mv \${outfile}.tmp \${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}

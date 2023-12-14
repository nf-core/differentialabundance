process FILTER_DIFFTABLE {

    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(input_file)
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
    output_file=\$(echo $input_file | sed 's/\\(.*\\)\\..*/\\1/')_filtered.tsv

    # Function to find column number
    find_column_number() {
        awk -v column="\$2" '{for(i=1;i<=NF;i++) if (\$i == column) {print i; exit}}' <<< "\$(head -n 1 "\$1")"
    }

    # Extract column numbers
    logFC_col=\$(find_column_number "$input_file" "log2FoldChange")
    padj_col=\$(find_column_number "$input_file" "padj")

    # Prepare the output file
    head -n 1 "$input_file" > "\${output_file}.tmp"

    # The following snippet performs the following checks on each row (add +0.0 to the numbers so that they are definitely treated as numerics):
    #
    # 1. Check that the current logFC/padj is not NA
    # 2. Check that the current logFC is >= threshold (abs does not work, so use a workaround)
    # 3. Check that the current padj is <= threshold
    #
    # If this is true, the row is written to the new file, otherwise not

    awk -F'\\t' -v logFC_col="\$logFC_col" -v padj_col="\$padj_col" -v logFC_thresh="$logFC_threshold" -v padj_thresh="$padj_threshold" '
        NR > 1 && \$logFC_col != "NA" && \$padj_col != "NA" &&
        ((\$logFC_col+0.0 >= logFC_thresh+0.0) || (-\$logFC_col+0.0 >= logFC_thresh+0.0)) &&
        \$padj_col+0.0 <= padj_thresh+0.0 { print }
    ' "$input_file" >> "\${output_file}.tmp"

    # Rename temporary file to final output file
    mv "\${output_file}.tmp" "\$output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}

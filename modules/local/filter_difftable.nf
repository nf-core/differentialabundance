process FILTER_DIFFTABLE {

    label 'process_single'

    conda "pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(input_file)
    tuple val(logFC_column), val(FC_threshold)
    tuple val(padj_column), val(padj_threshold)

    output:
    tuple val(meta), path("*_filtered.tsv") , emit: filtered
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    from math import log2
    from os import path
    import pandas as pd
    import platform
    from sys import exit

    # 1. Check that the current logFC/padj is not NA
    # 2. Check that the current logFC is >= threshold (abs does not work, so use a workaround)
    # 3. Check that the current padj is <= threshold
    # If this is true, the row is written to the new file, otherwise not
    if not any("$input_file".endswith(ext) for ext in [".csv", ".tsv", ".txt"]):
        exit("Please provide a .csv, .tsv or .txt file!")

    table = pd.read_csv("$input_file", sep=("," if "$input_file".endswith(".csv") else "\t"), header=0)
    logFC_threshold = log2(float("$FC_threshold"))
    table = table[~table["$logFC_column"].isna() &
                ~table["$padj_column"].isna() &
                (pd.to_numeric(table["$logFC_column"], errors='coerce').abs() >= float(logFC_threshold)) &
                (pd.to_numeric(table["$padj_column"], errors='coerce') <= float("$padj_threshold"))]

    table.to_csv(path.splitext(path.basename("$input_file"))[0]+"_filtered.tsv", sep="\t", index=False)

    with open('versions.yml', 'a') as version_file:
        version_file.write('"${task.process}":' + "\\n")
        version_file.write("    pandas: " + str(pd.__version__) + "\\n")
    """
}

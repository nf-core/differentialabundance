process TABULAR_TO_GSEA_CHIP {

    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path tsv
    tuple val(id), val(symbol)

    output:
    path "*.chip"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    function find_column_number {
        file=\$1
        column=\$2

        head -n 1 \$file | tr '\\t' '\\n' | grep -n "^\${column}\$" | awk -F':' '{print \$1}'
    }

    id_col=\$(find_column_number $tsv $id)
    symbol_col=\$(find_column_number $tsv $symbol)
    outfile=\$(echo $tsv | sed 's/\\(.*\\)\\..*/\\1/').chip

    echo -e "Probe Set ID\\tGene Symbol\\tGene Title" > \${outfile}.tmp
    tail -n +2 $tsv | awk -F'\\t' -v id=\$id_col -v symbol=\$symbol_col '{print \$id"\\t"\$symbol"\\tNA"}' >> \${outfile}.tmp
    mv \${outfile}.tmp \${outfile}
    """
}
